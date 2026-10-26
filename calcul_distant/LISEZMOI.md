# Calculs sur les runners GitHub (workflow `gazeification`)

Le workflow `.github/workflows/gazeification.yml` compile le solveur dans le conteneur `openfoam/openfoam8-paraview56`
(compilation mise en cache tant que `src/`, `applications/` et `Allwmake` ne changent pas), lance une matrice de cas,
arrête chaque passage proprement avant la limite de 6 h, sauvegarde les cas en artefacts et se relance lui-même
jusqu'à la fin des calculs.

## Lancer

```
gh workflow run gazeification.yml -R symphorientchimoe-ux/porousGasificationFoamMR \
   -f cas='["dossier/cas1", "dossier/cas2+dossier/cas3"]' -f ref=v1.0 -f np=2 -f etiquette=serie1
```

| Entrée | Rôle | Défaut |
|---|---|---|
| `cas` | liste JSON de dossiers du dépôt ; `a+b` fait tourner deux cas en même temps sur la même machine (4 vCPU : 2 processus chacun) | tutoriel `kimReactor149_restart` |
| `ref` | branche, étiquette ou commit des sources compilées | celui du workflow |
| `np` | processus MPI par cas (un cas déjà découpé garde son découpage) | 2 |
| `budget` | horloge [s] après laquelle le passage s'arrête au **prochain instant d'écriture** | 15000 |
| `budget_max` | horloge [s] d'arrêt immédiat avec écriture (filet de sécurité) | 19800 |
| `ecriture` | `writeInterval` imposé [s simulées] ; vide = celui du cas | vide |
| `reprise_run_id`, `passe` | gérés par le job « suite » | — |
| `passes_max` | nombre maximal de passages enchaînés | 40 |
| `etiquette` | suffixe des artefacts, pour distinguer des séries | vide |

Un fichier `ARRET` dans le cas (« colonne délai ») arrête le calcul « délai » secondes simulées après le passage
à 600 °C de la sonde Ts de cette colonne (`postProcessing/probesKim_fluid/fluidZone/*/Ts`, colonne 1 = temps).

## Règles

- OpenFOAM 8 compte les instants d'écriture et de sondes **depuis le temps de démarrage**. Le passage s'arrête
  donc au prochain instant d'écriture prévu (`nextWrite`), pour que chaque reprise parte d'un instant régulier.
  Choisir `writeInterval` multiple de l'intervalle des sondes, et assez court pour qu'un intervalle tienne dans
  `budget_max − budget` (avec le cas Kim à 5 mm, 40 à 120 s simulées par heure : 60 s simulées ≤ 1 h 30).
- Avec la v1.1 (`exactRestart yes`), une reprise suit le calcul continu (écart de 0,012 K sur T après 0,05 s,
  contre 8,5 K en v1.0) ; `lancer_cas.sh` écrit les temps à 17 chiffres pour cela.
- 20 machines au plus pour tout le compte GitHub (partagées avec les autres dépôts).

## Sorties

- `cas-<élément>` : archive du cas (dernier temps écrit) pour la reprise, gardée 7 jours ;
- `bilan-<élément>-passeN` : état (`termine`, `a_reprendre`, `erreur`), bilan (temps simulé, vitesse, passages
  des thermocouples), fin du journal, gardé 90 jours. Résumé aussi sur la page du run.

## Coût mesuré (26/10/2026, tutoriel en série)

Compilation 3 min (puis cache). Passage : 2 à 2,5 min de surcoût (conteneur, artefacts). Vitesse en série :
31 à 50 s simulées par heure selon la machine (station : 42 s/h sur un cœur libre).
