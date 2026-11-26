FoamFile { version 2.0; format ascii; class dictionary; object topoSetDict; }
// 2026-10-06 Kim : couvercle du reacteur (Fig. 1 : reacteur ferme, fumees par une conduite centrale vers
// condenseurs + extraction). Faces de sortie r > 0.02 m -> paroi "lid" ; ouverture centrale d = 40 mm
// (diametre NON donne par Kim : HYPOTHESE declaree).
actions
(
    { name lidFaces; type faceSet; action new;    source patchToFace; sourceInfo { name outlet; } }
    { name lidFaces; type faceSet; action subset; source boxToFace; sourceInfo { box (0.02 -1 0.69) (1 1 0.71); } }
);
