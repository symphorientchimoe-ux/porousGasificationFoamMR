#!/usr/bin/env python3
"""Campagne v1.1 (rapport d'expertise CPC, étape B, 26/10/2026) : construit les cas de la vague 1 à partir des états
à 60 s de l'étape 1 (base_<débit>, v1.1 + cinétique de Mehrabian), dans DEST/<nom>.
Chaque cas = base_<débit> (constant, system, processor*/constant, processor*/60) + modifications listées ci-dessous.
Étape 2 : allumeur 1000 K de 60 à 300 s (sauf indication), maxUa 0,3 (comme avant), writePrecision 17.
Couvercle : maillage de case_woodPellet_149_Kim_ref5_lid, champs à 60 s adaptés par _src60/poser_couvercle.py
(même méthode que le 16/10, CONSTATS 21.41).
Usage : python3 construire.py DEST [nom ...]"""
import os, re, shutil, subprocess, sys

ICI = os.path.dirname(os.path.abspath(__file__))
LID = os.path.join(ICI, '..', 'case_woodPellet_149_Kim_ref5_lid')
POSER = os.path.join(LID, '_src60', 'poser_couvercle.py')
BUILD = '/home/cndt/OpenFOAM/work/porousGasificationFoam_OF8/release/_build_v11b'

def sens(fichier, motif, facteur):
    """multiplie par facteur le nombre qui suit le motif (première occurrence hors commentaire)"""
    return (fichier, motif, facteur)

# nom : (débit, couvercle, modifications, ARRET (colonne, délai) ou None, endTime)
# colonnes des sondes Ts : 2 = T0 (z 0,1005), 3 = T1, 4 = T2, 5 = T3, 6 = T4, 7 = T5
HV = 'constant/fluidZone/heatTransferProperties'; HM = 'constant/fluidZone/specieTransferProperties'
RAD = 'constant/fluidZone/radiationProperties'; KR = 'constant/solidZone/thermophysicalProperties'
STP = 'constant/fluidZone/solidThermophysicalProperties'; CH = 'constant/fluidZone/chemistryProperties'
CB = 'constant/fluidZone/combustionProperties'; PY = 'constant/fluidZone/pyrolysisProperties'
CAS = {
    # références et essais en aveugle (avec couvercle, comme l'expérience)
    'R149L':  (149, True,  [], (2, 420), 9000),          # référence, jusqu'à T0 + 420 s (fin de propagation)
    'R149O':  (149, False, [], None, 15000),             # ouvert, essai complet (250 min de Kim) : gaz de sortie, phase charbon
    'L74':    (74,  True,  [], (3, 420), 9000),
    'L298':   (298, True,  [], (3, 420), 6000),
    'L596':   (596, True,  [], (3, 420), 6000),
    'O596':   (596, False, [], (3, 420), 4000),          # diagnostic : ouvert (comme le 596 v1.0 qui ne s'allumait pas)
    'L596i600': (596, True, [(PY, r'igniterDuration\s+', 600 / 300)], (3, 420), 6000),   # allumeur tenu 600 s
    # ensemble (incertitude numérique) : débit d'air perturbé de 1e-8 en relatif
    'E1':     (149, False, [('U_inlet', None, 1 + 1e-8)], (3, 420), 6000),
    'E2':     (149, False, [('U_inlet', None, 1 - 1e-8)], (3, 420), 6000),
    'E3':     (149, False, [('U_inlet', None, 1 + 2e-8)], (3, 420), 6000),
    # sensibilités (ouvert, base = R149O), arrêt 60 s après T2 (trois intervalles T5 -> T2)
    'S_hv05': (149, False, [(HV, r'\n\s*h\s+', 0.5), (HM, r'\n\s*h\s+', 0.5)], (4, 60), 6000),
    'S_hv2':  (149, False, [(HV, r'\n\s*h\s+', 2.0), (HM, r'\n\s*h\s+', 2.0)], (4, 60), 6000),
    'S_as05': (149, False, [(RAD, r'\n\s*as\s+as\s+\[[^\]]*\]\s+', 0.5), (RAD, r'\n\s*borderAs\s+borderAs\s+\[[^\]]*\]\s+', 0.5)], (4, 60), 6000),
    'S_as2':  (149, False, [(RAD, r'\n\s*as\s+as\s+\[[^\]]*\]\s+', 2.0), (RAD, r'\n\s*borderAs\s+borderAs\s+\[[^\]]*\]\s+', 2.0)], (4, 60), 6000),
    'S_kr05': (149, False, [(KR, r'\n\s*kappa\s+', 0.5)], (4, 60), 6000),
    'S_kr2':  (149, False, [(KR, r'\n\s*kappa\s+', 2.0)], (4, 60), 6000),
    'S_rc244': (149, False, [(STP, r'density\s+\{\s*rho\s+(?=305)', 244 / 305)], (4, 60), 6000),
    'S_rc488': (149, False, [(STP, r'density\s+\{\s*rho\s+(?=305)', 488 / 305)], (4, 60), 6000),
    'S_dvm':  (149, False, [(CH, r'\n\s*\((?=4\.4e9 16982)', 0.0550232)], (4, 60), 6000),   # x exp(-2,9)
    'S_dvp':  (149, False, [(CH, r'\n\s*\((?=4\.4e9 16982)', 18.1741454)], (4, 60), 6000),  # x exp(+2,9)
    'S_edm05': (149, False, [(CB, r'\n\s*A\s+', 0.5 / 4)], (4, 60), 6000),
    'S_edmInf': (149, False, [(CB, r'\n\s*A\s+', 1e6 / 4)], (4, 60), 6000),             # A = 1e6 : chimie cinétique seule
    'S_sa06': (149, False, [(PY, r'\n\s*shrinkA\s+', 0.6 / 0.8)], (4, 60), 6000),
    'S_sa10': (149, False, [(PY, r'\n\s*shrinkA\s+', 1.0 / 0.8)], (4, 60), 6000),
    # maillage (GCI, point M8) : 2,5 et 10 mm, état v1.1 à 60 s transféré par mapFields (mapNearest), arrêt 60 s après T3
    'M25':    ('149_M25', False, [], (5, 60), 6000),
    'M10':    ('149_M10', False, [], (5, 60), 6000),
}

def multiplier(path, motif, facteur):
    # newline='' : le fichier U à 60 s est en binaire, la conversion des fins de ligne le corrompait (E1-E3, 26/10)
    s = open(path, encoding='latin1', newline='').read()
    m = re.search(motif + r'([0-9.eE+-]+)', s)
    if not m:
        sys.exit(f'motif introuvable : {motif} dans {path}')
    v = float(m.group(1)) * facteur
    s = s[:m.start(1)] + f'{v:.10g}' + s[m.end(1):]
    open(path, 'w', encoding='latin1', newline='').write(s)
    return m.group(1), f'{v:.10g}'

def regler(path, cle, valeur):
    s = open(path, encoding='latin1').read()
    s, n = re.subn(rf'(?m)^({cle}\s+)[^;]+;.*$', rf'\g<1>{valeur};', s)
    if n == 0:
        s += f'\n{cle} {valeur};\n'
    open(path, 'w', encoding='latin1').write(s)

def construire(dest, nom):
    deb, lid, mods, arret, fin = CAS[nom]
    base = os.path.join(ICI, f'base_{deb}')
    d = os.path.join(dest, nom)
    shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    for x in ('constant', 'system'):
        shutil.copytree(os.path.join(base, x), os.path.join(d, x), ignore=shutil.ignore_patterns('*.bak*', 'JOURNAL_PROMPTS.md'))
    procs = sorted(p for p in os.listdir(base) if p.startswith('processor'))
    for p in procs:
        shutil.copytree(os.path.join(base, p, 'constant'), os.path.join(d, p, 'constant'))
        shutil.copytree(os.path.join(base, p, '60'), os.path.join(d, p, '60'))
    journal = [f'{nom} : débit {deb} kg/m2h, couvercle {"oui" if lid else "non"}, départ base_{deb} à 60 s']
    if lid:
        # même suite que le 16/10 (CONSTATS 21.41) : reconstruction à 60 s, champs en ascii, maillage avec couvercle,
        # poser_couvercle.py (sortie tronquée à 2 faces + bloc lid repris de ref5_lid), puis redécoupage
        cmd = f"""source ~/OpenFOAM/OpenFOAM-8/etc/bashrc > /dev/null 2>&1
set -e
export LD_LIBRARY_PATH={BUILD}/lib:$LD_LIBRARY_PATH
cd {d}
reconstructPar -allRegions -time 60 > log.reconstruct 2>&1
sed -i -E 's|^(writeFormat\\s+)\\w+;|\\1ascii;|' system/controlDict
foamFormatConvert -region fluidZone -time 60 > log.ffc 2>&1
rm -rf constant/fluidZone/polyMesh; cp -r {LID}/constant/fluidZone/polyMesh constant/fluidZone/
for x in 60/fluidZone/*; do [ -e {LID}/_src60/60/fluidZone/$(basename $x) ] || [ $(basename $x) = phi ] || rm -f $x; done
python3 {POSER} 60/fluidZone {LID}/_src60/60/fluidZone > log.couvercle 2>&1
rm -f 60/fluidZone/phi.retire
rm -rf processor*
decomposePar -allRegions -time 60 > log.decompose 2>&1
rm -rf 60
"""
        r = subprocess.run(['bash', '-c', cmd], capture_output=True, text=True)
        if r.returncode:
            sys.exit(f'couvercle {nom} : {r.stdout} {r.stderr}')
        procs = sorted(p for p in os.listdir(d) if p.startswith('processor'))
        journal.append('couvercle posé à 60 s (poser_couvercle.py)')
    for f, motif, fac in mods:
        if f == 'U_inlet':
            for p in procs:
                path = os.path.join(d, p, '60', 'fluidZone', 'U')
                s = open(path, encoding='latin1', newline='').read()
                if 'massFlowRate' not in s:
                    continue
                a, b = multiplier(path, r'massFlowRate\s+constant\s+', fac)
            journal.append(f'débit d air x {fac:.10g} ({a} -> {b})')
        else:
            a, b = multiplier(os.path.join(d, f), motif, fac)
            journal.append(f'{f} : {a} -> {b}')
    py = os.path.join(d, PY)
    s = open(py, encoding='latin1').read()
    s = re.sub(r'(?m)^(    igniterTemperature  )[0-9]+;', r'\g<1>1000;', s)
    open(py, 'w', encoding='latin1').write(s)
    cd = os.path.join(d, 'system', 'controlDict')
    regler(cd, 'endTime', fin); regler(cd, 'maxUa', 0.3); regler(cd, 'writePrecision', 17)
    regler(cd, 'startFrom', 'latestTime'); regler(cd, 'stopAt', 'endTime')
    if arret:
        open(os.path.join(d, 'ARRET'), 'w').write(f'{arret[0]} {arret[1]}\n')
        journal.append(f'arrêt {arret[1]} s après le passage à 600 °C de la colonne {arret[0]}')
    journal.append(f'endTime {fin}')
    open(os.path.join(d, 'LISEZMOI.txt'), 'w').write('\n'.join(journal) + '\n')
    print(' | '.join(journal))

if __name__ == '__main__':
    dest = sys.argv[1]
    for nom in (sys.argv[2:] or CAS):
        construire(dest, nom)
