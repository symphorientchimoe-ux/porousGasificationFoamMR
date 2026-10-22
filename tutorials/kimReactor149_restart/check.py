#!/usr/bin/env python3
"""Compares a run of kimReactor149_restart with reference/: last probe values (relative difference), bed-top heights
after each collapse, and the conservation lines of the collapse remap. Prints PASS or FAIL for each check.
Usage: python3 check.py [log file]   (default log.porousGasificationFoamMR)"""
import sys, os, re, glob
import numpy as np

log = sys.argv[1] if len(sys.argv) > 1 else 'log.porousGasificationFoamMR'
ok = True

def last_row(path):
    rows = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    rows = [r for r in rows if not any('(' in x for x in r)]   # scalar probes only
    return np.array(rows[-1], float) if rows else None

# 1. probes: last row (t = 2402 s), relative difference
worst = (0.0, '')
ncmp = 0
for ref in sorted(glob.glob('reference/postProcessing/probes*/**/*', recursive=True)):
    if os.path.isdir(ref):
        continue
    new = ref[len('reference/'):]
    if not os.path.exists(new):
        print(f'FAIL  missing output {new}'); ok = False; continue
    a, b = last_row(ref), last_row(new)
    if a is None:
        continue
    if b is None or a.shape != b.shape or abs(a[0] - b[0]) > 1e-6:
        print(f'FAIL  {new}: last time or size differs'); ok = False; continue
    ncmp += 1
    d = np.max(np.abs(b[1:] - a[1:]) / np.maximum(np.abs(a[1:]), 1e-12))
    if d > worst[0]:
        worst = (d, new)
tol = 1e-3
print(f'{"PASS" if worst[0] < tol and ncmp else "FAIL"}  probes at t = 2402 s ({ncmp} files): max relative difference {worst[0]:.2e}, tolerance {tol:g}')
ok &= worst[0] < tol and ncmp > 0

# 2. bed collapse: heights and conservation
txt = open(log).read()
hr = [float(x) for x in re.findall(r'bed top axis column : \S+ -> (\S+) m', open('reference/bedCollapse.txt').read())]
hn = [float(x) for x in re.findall(r'bed top axis column : \S+ -> (\S+) m', txt)]
if len(hn) != len(hr):
    print(f'FAIL  number of collapses {len(hn)} (reference {len(hr)})'); ok = False
else:
    dh = max(abs(x - y) for x, y in zip(hn, hr))
    print(f'{"PASS" if dh < 1e-6 else "FAIL"}  bed-top height on the axis after {len(hn)} collapses: max difference {dh:.1e} m')
    ok &= dh < 1e-6
rel = [abs(float(x)) for x in re.findall(r'conservation .* rel (\S+)', txt)]
mr = max(rel) if rel else float('nan')
print(f'{"PASS" if rel and mr < 1e-12 else "FAIL"}  mass conservation of the remap: max relative error {mr:.1e} ({len(rel)} lines)')
ok &= bool(rel) and mr < 1e-12
print('ALL CHECKS PASSED' if ok else 'SOME CHECKS FAILED')
sys.exit(0 if ok else 1)
