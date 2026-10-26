#!/usr/bin/env python3
"""Court bilan d'un cas porousGasificationFoamMR : temps simulé, vitesse du dernier passage,
passages à 600 °C des sondes Ts (postProcessing/probesKim_fluid/fluidZone/*/Ts), état.
Usage : bilan_cas.py DOSSIER_CAS [--jalon COLONNE DELAI]
Avec --jalon, écrit DOSSIER_CAS/ARRET_ATTEINT quand le temps simulé dépasse le passage à 600 °C
de la colonne COLONNE (1 = temps) de DELAI secondes."""
import glob, os, re, sys

case = sys.argv[1]
jalon = None
if '--jalon' in sys.argv:
    i = sys.argv.index('--jalon')
    jalon = (int(sys.argv[i + 1]), float(sys.argv[i + 2]))
SEUIL = 873.15

def logs():
    return sorted(glob.glob(os.path.join(case, 'log.run.*')), key=lambda p: int(p.rsplit('.', 1)[1]))

def temps_log(p):
    t = re.findall(r'^Time = (\S+)', open(p, errors='replace').read(), re.M)
    return (float(t[0]), float(t[-1]), len(t)) if t else None

def sondes():
    rows = {}
    files = glob.glob(os.path.join(case, 'postProcessing/probesKim_fluid/fluidZone/*/Ts'))
    for f in sorted(files, key=lambda p: float(p.split('/')[-2])):
        for l in open(f):
            if l.startswith('#'):
                continue
            v = l.split()
            try:
                rows[float(v[0])] = [float(x) for x in v[1:]]
            except ValueError:
                pass
    return [(t, rows[t]) for t in sorted(rows)]

def passage(rows, col):
    i, prev = col - 2, None
    for t, r in rows:
        if i < len(r) and r[i] > SEUIL:
            if prev is None:
                return t
            t0, v0 = prev
            return t0 + (SEUIL - v0) * (t - t0) / (r[i] - v0)
        prev = (t, r[i]) if i < len(r) else None
    return None

rows = sondes()
if jalon:
    tj = passage(rows, jalon[0])
    if tj is not None and rows and rows[-1][0] >= tj + jalon[1]:
        open(os.path.join(case, 'ARRET_ATTEINT'), 'w').write(f'{tj:.1f} {jalon[1]:g}\n')
    sys.exit(0)

L = logs()
if L:
    tl = temps_log(L[-1])
    duree = float(open(os.path.join(case, 'duree_passage.txt')).read()) if os.path.exists(os.path.join(case, 'duree_passage.txt')) else None
    if tl:
        print(f'passage {len(L)} : t = {tl[0]:.3f} -> {tl[1]:.3f} s, {tl[2]} pas'
              + (f', {duree:.0f} s d horloge, {3600 * (tl[1] - tl[0]) / max(duree, 1):.1f} s simulees par heure' if duree else ''))
if rows:
    noms = {7: 'T5', 6: 'T4', 5: 'T3', 4: 'T2', 3: 'T1'}
    print('passages a 600 C (sondes Ts) : ' + ', '.join(
        f'{n} ' + (f'{passage(rows, c):.1f} s' if passage(rows, c) is not None else '-') for c, n in sorted(noms.items(), reverse=True)))
etat = os.path.join(case, 'ETAT')
if os.path.exists(etat):
    print('etat : ' + open(etat).read().strip())
