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
