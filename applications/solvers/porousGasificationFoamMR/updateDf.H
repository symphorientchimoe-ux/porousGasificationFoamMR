/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  | Modified work Copyright (C) 2026 S. Kemle Tchimoe et al., CNDT/MINRESI
-------------------------------------------------------------------------------
License
    This file is part of porousGasificationFoamMR, a derivative of OpenFOAM and
    of porousGasificationFoam (Zuk et al., Comput. Phys. Commun. 278 (2022) 108407).

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// Darcy resistance Df (Ergun) recomputed at every time step from the local,
// evolving porosity, in the cellZone biomass_bed only. In porousGasificationFoam
// Df is read once at start-up (createFields.H) and never updated, while porosityF
// (the same field as the porosity_ of volPyrolysis, see createPyrolysisModel.H)
// rises towards 1 as a cell is consumed: without this update the gas keeps the
// resistance of the initial dense bed even in a fully burnt cell.
//
// Same formula as the initial Df of system/setFieldsDict:
//   Df(eps) = 150*(1-eps)^2 / (eps^3 * dp^2)
// Continuous at t = 0: eps = 0.4, dp = 0.05 m gives exactly 3.375e5.
//
// eps is clamped to [0.05, 0.99] (a numerical choice): avoids Df = 0 (no
// resistance, possible local velocity peak) when porosityF -> 1, and the
// singularity at eps -> 0.
//
// Cells outside biomass_bed keep the value set in system/setFieldsDict.
{
    const label biomassBedZoneID = mesh.cellZones().findZoneID("biomass_bed");

    if (biomassBedZoneID != -1)
    {
        const labelList& bedCells = mesh.cellZones()[biomassBedZoneID];

        // Opt-in (2026-09-27) : Ergun inertial term evaluated with the LOCAL
        // porosity, folded into Df as 1.75(1-eps)/(eps^3 dp)*rho|U|/mu (lagged U),
        // instead of a constant Forchheimer coefficient computed for eps=0.481
        // (constant/fluidZone/porosityProperties must then set forchheimerCoeff 0).
        // controlDict: ergunInertialEps yes; ergunDp 0.049;  (defaults: no, 0.05)
        const bool ergunInertialEps =
            runTime.controlDict().lookupOrDefault<Switch>("ergunInertialEps", false);
        const scalar dp =
            runTime.controlDict().lookupOrDefault<scalar>("ergunDp", 0.05);
        tmp<volScalarField> tmuErgun(thermo.mu());
        const volScalarField& muErgun = tmuErgun();

        forAll(bedCells, i)
        {
            const label celli = bedCells[i];

            const scalar eps = max(min(porosityF[celli], 0.99), 0.05);

            scalar DfVal = 150.0*sqr(1.0 - eps)/(pow3(eps)*sqr(dp));
            if (ergunInertialEps)
            {
                DfVal +=
                    1.75*(1.0 - eps)/(pow3(eps)*dp)
                   *rho[celli]*mag(U[celli])/max(muErgun[celli], small);
            }

            Df[celli] = tensor
            (
                DfVal, 0, 0,
                0, DfVal, 0,
                0, 0, DfVal
            );
        }

        Df.correctBoundaryConditions();
    }
}
