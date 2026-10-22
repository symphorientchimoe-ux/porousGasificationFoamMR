/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object topoSetDict; }

actions
(
    { name lidFaces; type faceSet; action new;    source patchToFace; sourceInfo { name outlet; } }
    { name lidFaces; type faceSet; action subset; source boxToFace; sourceInfo { box (0.02 -1 0.69) (1 1 0.71); } }
);

// ************************************************************************* //
