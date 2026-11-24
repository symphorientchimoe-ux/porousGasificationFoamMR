# Série de maillage avec mélange de pore (21/11/2026)

Article ECM, relectures 1 et 2 : la vitesse du front ne converge pas en maillage avec le mélange de l'EDM fondé sur le
taux de déformation résolu (sqrt(2)|symm(grad U)|), qui croît quand on raffine. Option `mixing poreScale` de laminarEDM :
dans le lit, fréquence de traversée d'un pore |U|/(eps dp) + (mu/rho)/dp^2, sans gradient. Défaut inchangé (v1.1 au bit près,
vérifié sur kimReactor149_restart).

| Cas | maille axiale dans le lit | arrêt | processus |
|---|---|---|---|
| P10 | 10 mm | T3 + 60 s | 2 |
| P5 | 5 mm | T3 + 60 s | 2 |
| P25 | 2,5 mm | T3 + 60 s | 2 |
| P125 | 1,25 mm | T4 + 60 s | 4 |

149 kg/m2h, ciel ouvert, allumeur 10 mm, état à 60 s identique aux séries M (mapFields depuis 5 mm).
| Q5 | 5 mm, effondrement coupé | T4 + 60 s | 2 |
| Q25 | 2,5 mm, effondrement coupé | T4 + 60 s | 2 |

## Cas principaux de l'article ECM sur maille de 2,5 mm (24/11/2026, option 2)
P5 / P25 / P125 : T5 -> T4 = 862,9 / 805,3 / 800,5 s, convergé dès 2,5 mm avec le mélange de pore. Les cas principaux
(avec couvercle) sont donc refaits à 2,5 mm. Construction : validation_Kim2016/campagne_v11/maillage25/construire_L25.sh
(maillage et réglages de P25, patch lid comme ref5_lid, état à 60 s des cas 5 mm transféré par mapFields mapNearest,
propriétés de transfert et débit de la référence). Essai local de L25_149 : démarrage normal (60 -> 62 s).

| Cas | débit (kg/m2h) | référence 5 mm (état à 60 s) | arrêt |
|---|---|---|---|
| L25_149 | 149 | R149L (+ bilanPostes et coupes, comme B149bilan) | T1 + 420 s |
| L25_74 | 74,5 | L74 | T1 + 420 s |
| L25_298 | 298 | L298 | T1 + 420 s |
| L25_596 | 596 | L596_ign4r (allumeur 40 mm, igniterZone refait ; 60-65 s en local) | T1 + 420 s |
