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

scalar DiNum =  pyrolysis->solidRegionDiffNo();

// 2026-09-09: adds the diffusion number of the solidZone region (refractory,
// pure conduction) to the global time-step control, see createSolidFields.H
DiNum = max
(
    DiNum,
    solidRegionDiffNo
    (
        solidMesh,
        runTime,
        rhoSolid*thermoSolid.Cp()(),
        thermoSolid.kappa()()
    )
);

// 2026-09-17: radiative relaxation number,
// same lagged-estimate principle as CoNum/DiNum (uses radiationRp_/T_ as
// they stood after the PREVIOUS radiation.H call, this control runs
// before radiation.H/evolve() each timestep). Folded into
// setMultiRegionDeltaT.H via its own maxRa (separate from maxDi -- the
// radiative sink can overshoot even under implicit treatment when the
// linearisation goes stale over one large step, so it needs a stricter
// cap than the diffusion number).
scalar RaNum = pyrolysis->radiativeRelaxNo();

// 2026-09-21 : gas-phase relaxation number, same lagged-estimate
// principle as RaNum above but for the explicit gas-phase energy source
// terms (heatTransfer*porosityF + heatUpGas) in EEqn.H. Folded into
// setMultiRegionDeltaT.H via its own maxGa.
scalar GaNum = pyrolysis->gasRelaxNo();

// 2026-09-21 : velocity-growth number, same lagged-estimate principle as
// CoNum/DiNum/RaNum/GaNum but for the momentum equation directly --
// U at this point is the value at the end of the PREVIOUS timestep,
// U.oldTime() the one before that (Euler ddt already keeps it, no new
// stored field needed), so this measures how much U changed, as a
// fraction of itself, over the last dt -- exactly the growth that
// Courant (U*dt/dx, a spatial-transport measure) does not catch. Folded
// into setMultiRegionDeltaT.H via its own maxUa. See
// readPyrolysisTimeControls.H for the full rationale.
// 2026-09-29: velocity floor of the Ua criterion read from controlDict
// (UaSmall, default 1e-3 m/s = previous behaviour). In the Kim case the air below the
// grate moves at 0.2-2 mm/s and its relative acoustic oscillations exceed 0.3, giving dt ~1e-5 s.
scalar UaNum = gMax
(
    (
        mag(U - U.oldTime())
      / (mag(U.oldTime()) + dimensionedScalar("smallU", dimVelocity, runTime.controlDict().lookupOrDefault<scalar>("UaSmall", 1e-3)))
    )()
);
