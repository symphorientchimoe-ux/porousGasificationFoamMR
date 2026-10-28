#!/bin/bash
# Lance (ou reprend) un cas porousGasificationFoamMR (OpenFOAM 8) sur n'importe quelle machine
# (poste local, conteneur Docker d'un runner GitHub).
#
# Le calcul s'arrête proprement (écriture des champs) au bout de BUDGET secondes d'horloge.
# En relançant la même commande sur le même dossier, il repart du dernier temps écrit.
#
# Usage :
#   lancer_cas.sh DOSSIER_CAS [--np 2] [--budget 15000] [--budget-max 19800] [--ecriture DT]
#                 [--solveur porousGasificationFoamMR] [--archive sortie.tar.gz]
#
# Au bout de --budget s d'horloge, le calcul s'arrête au PROCHAIN instant d'écriture prévu (nextWrite) :
# OpenFOAM 8 compte les instants d'écriture (champs et sondes) depuis le temps de démarrage, si bien
# qu'un arrêt à un instant quelconque décalerait toutes les écritures des passages suivants (et
# endTime ne serait plus écrit). Filet de sécurité : arrêt immédiat avec écriture à --budget-max s.
# --ecriture DT fixe writeInterval (s simulées) : à choisir pour qu'un intervalle tienne dans
# budget-max - budget.
#
# Arrêt sur jalon (facultatif) : un fichier DOSSIER_CAS/ARRET contenant « colonne délai »
# arrête le calcul « délai » secondes simulées après le premier passage à 873,15 K (600 °C)
# de la sonde Ts de cette colonne (colonne 1 = temps) dans postProcessing/probesKim_fluid/fluidZone/*/Ts.
#
# À la fin, DOSSIER_CAS/ETAT contient « termine », « a_reprendre » ou « erreur »,
# et DOSSIER_CAS/bilan.txt un court bilan (temps simulé, vitesse, passages des thermocouples).
#
# Variables d'environnement : FOAM_BASHRC (défaut : /opt/openfoam8/etc/bashrc s'il existe,
# sinon ~/OpenFOAM/OpenFOAM-8/etc/bashrc) ; PGF_BUILD (dossier contenant bin/ et lib/ du solveur compilé).

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -z "${FOAM_BASHRC:-}" ]; then
    FOAM_BASHRC=$(ls -1 /opt/openfoam8/etc/bashrc 2>/dev/null | head -1)
    [ -n "$FOAM_BASHRC" ] || FOAM_BASHRC="$HOME/OpenFOAM/OpenFOAM-8/etc/bashrc"
fi

CASE="$1"; shift
NP=2; BUDGET=15000; BUDGET_MAX=19800; ECRITURE=""; SOLVEUR=porousGasificationFoamMR; ARCHIVE=""
while [ $# -gt 0 ]; do
    case "$1" in
        --np) NP="$2"; shift 2 ;;
        --budget) BUDGET="$2"; shift 2 ;;
        --budget-max) BUDGET_MAX="$2"; shift 2 ;;
        --ecriture) ECRITURE="$2"; shift 2 ;;
        --solveur) SOLVEUR="$2"; shift 2 ;;
        --archive) ARCHIVE="$2"; shift 2 ;;
        *) echo "option inconnue $1"; exit 2 ;;
    esac
done
CASE="$(realpath -m "$CASE")"
[ -n "$ARCHIVE" ] && ARCHIVE="$(realpath -m "$ARCHIVE")"

# OpenFOAM refuse mpirun en root sans ces variables (conteneurs)
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export OMPI_MCA_btl_vader_single_copy_mechanism=none

set +u; source "$FOAM_BASHRC"; set -u
if [ -n "${PGF_BUILD:-}" ]; then
    export PATH="$PGF_BUILD/bin:$PATH" LD_LIBRARY_PATH="$PGF_BUILD/lib:${LD_LIBRARY_PATH:-}"
fi
command -v "$SOLVEUR" > /dev/null || { echo "solveur $SOLVEUR introuvable"; echo erreur > "$CASE/ETAT"; exit 1; }
# comme sur la station : pas de piège à exceptions ni d'initialisation à NaN
unset FOAM_SIGFPE FOAM_SETNAN

cd "$CASE" || exit 1
rm -f ETAT ARRET_ATTEINT

# découpage si besoin
NPROC=$(ls -d processor[0-9]* 2>/dev/null | wc -l)
if [ "$NP" -gt 1 ] && [ "$NPROC" -eq 0 ]; then
    foamDictionary system/decomposeParDict -entry numberOfSubdomains -set "$NP" > /dev/null
    for r in system/*/; do
        [ -f "$r/decomposeParDict" ] && foamDictionary "$r/decomposeParDict" -entry numberOfSubdomains -set "$NP" > /dev/null
    done
    decomposePar -allRegions > log.decomposePar 2>&1 || { echo erreur > ETAT; exit 1; }
    NPROC=$NP
fi
[ "$NPROC" -gt 0 ] && NP=$NPROC

# reprise exacte : temps et pas de temps écrits à 17 chiffres (uniform/time), noms de dossiers lisibles à la reprise
foamDictionary system/controlDict -entry writePrecision -set 17 > /dev/null
foamDictionary system/controlDict -entry timePrecision -set 12 > /dev/null
foamDictionary system/controlDict -entry stopAt -set endTime > /dev/null
[ -n "$ECRITURE" ] && foamDictionary system/controlDict -entry writeInterval -set "$ECRITURE" > /dev/null
# arrêt propre : au prochain instant d'écriture après BUDGET s d'horloge, au plus tard à BUDGET_MAX s
for fo in stopAtClock stopAtClockMax; do
    foamDictionary system/controlDict -entry functions/$fo -remove > /dev/null 2>&1
done
foamDictionary system/controlDict -entry functions/stopAtClock -add \
    "{ type stopAtClockTime; libs (\"libutilityFunctionObjects.so\"); stopTime $BUDGET; action nextWrite; }" > /dev/null
foamDictionary system/controlDict -entry functions/stopAtClockMax -add \
    "{ type stopAtClockTime; libs (\"libutilityFunctionObjects.so\"); stopTime $BUDGET_MAX; action writeNow; }" > /dev/null

n=1; while [ -e "log.run.$n" ]; do n=$((n+1)); done
LOG="log.run.$n"

# surveillance du jalon (fichier ARRET) : demande writeNow au solveur (runTimeModifiable)
surveiller() {
    read -r COL DELAI < ARRET
    while sleep 60; do
        python3 "$HERE/bilan_cas.py" . --jalon "$COL" "$DELAI" > /dev/null 2>&1
        if [ -f ARRET_ATTEINT ]; then
            foamDictionary system/controlDict -entry stopAt -set writeNow > /dev/null
            return
        fi
    done
}
WATCH=""
if [ -f ARRET ]; then surveiller & WATCH=$!; fi

echo ">> $SOLVEUR sur $NP processus, budget ${BUDGET} s, journal $LOG"
t0=$(date +%s)
if [ "$NP" -gt 1 ]; then
    # --oversubscribe : les runners annoncent des vCPU (hyperthreads) que MPI ne compte pas comme cœurs
    mpirun --oversubscribe --bind-to none -np "$NP" "$SOLVEUR" -parallel 2>&1 | awk -v N=20 -f "$HERE/logfilter.awk" > "$LOG"
    rc=${PIPESTATUS[0]}
else
    "$SOLVEUR" 2>&1 | awk -v N=20 -f "$HERE/logfilter.awk" > "$LOG"
    rc=${PIPESTATUS[0]}
fi
# 2026-10-28 : journal filtré (un pas sur 20 + toutes les alertes, logfilter.awk comme sur la station).
# Le journal brut (~4 ko par pas, ~1 Go par passage) s'accumulait dans l'archive de reprise : 3-4 Go par cas
# à la passe 6, et l'archive S_hv tronquée à la passe 5 (« Unexpected EOF »).
[ -n "$WATCH" ] && kill "$WATCH" 2> /dev/null
echo "$(( $(date +%s) - t0 ))" > duree_passage.txt
echo ">> $SOLVEUR terminé (code $rc) en $(cat duree_passage.txt) s"
tail -3 "$LOG"

if [ "$NP" -gt 1 ]; then TDIR=processor0; else TDIR=.; fi
LAST=$(ls -1 "$TDIR" | grep -E '^[0-9.eE+-]+$' | sort -g | tail -1)
END=$(foamDictionary system/controlDict -entry endTime -value | sed 's/[;[:space:]].*//')
if [ "$rc" -ne 0 ]; then
    echo erreur > ETAT
elif [ -f ARRET_ATTEINT ]; then
    echo termine > ETAT
# fin atteinte : dernier temps écrit ou dernier pas calculé (après une reprise, OpenFOAM 8 compte les
# instants d'écriture depuis le temps de reprise, si bien que endTime n'est pas toujours écrit)
elif LASTRUN=$(grep -E '^Time = ' "$LOG" | tail -1 | awk '{print $3}'); \
     python3 -c "import sys; sys.exit(0 if max(float('$LAST'), float('${LASTRUN:-0}')) >= float('$END')*(1 - 1e-9) else 1)"; then
    echo termine > ETAT
else
    echo a_reprendre > ETAT
fi
echo ">> dernier temps écrit : $LAST / $END -> $(cat ETAT)"

python3 "$HERE/bilan_cas.py" . > bilan.txt 2>&1 || true
cat bilan.txt

# archive légère pour la reprise : seul le dernier temps écrit est gardé (0/ à la racine aussi)
if [ -n "$ARCHIVE" ]; then
    shopt -s nullglob
    for d in processor*/ .; do
        for t in $(ls -1 "$d" | grep -E '^[0-9.eE+-]+$' | sort -g | head -n -1); do
            [ "$t" = "0" ] && [ "$d" = "." ] && continue
            rm -rf "${d%/}/$t"
        done
    done
    cd ..
    tar czf "$ARCHIVE.tmp" "$(basename "$CASE")" && mv -f "$ARCHIVE.tmp" "$ARCHIVE"
    echo ">> archive : $ARCHIVE ($(du -h "$ARCHIVE" | cut -f1))"
fi
exit 0
