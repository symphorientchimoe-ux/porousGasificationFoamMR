/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Global
    setMultiRegionDeltaT

Description
    Reset the timestep to maintain a constant maximum courant and
    diffusion Numbers. Reduction of time-step is immediate, but
    increase is damped to avoid unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    if (CoNum == -great)
    {
        CoNum = small;
    }

    if (DiNum == -great)
    {
        DiNum = small;
    }

    if (RaNum == -great)
    {
        RaNum = small;
    }

    if (GaNum == -great)
    {
        GaNum = small;
    }

    if (UaNum == -great)
    {
        UaNum = small;
    }

    scalar maxDeltaTFluid = maxCo/(CoNum + small);
    scalar maxDeltaTSolid = maxDi/(DiNum + small);
    // 2026-09-17: radiative relaxation
    // constraint, same damped-growth treatment as the other two.
    scalar maxDeltaTRad = maxRa/(RaNum + small);
    // 2026-09-21 : gas-phase relaxation constraint (see gasRelaxNo() in
    // heterogeneousPyrolysisModel.H), same damped-growth treatment.
    scalar maxDeltaTGas = maxGa/(GaNum + small);
    // 2026-09-21 : velocity-growth constraint (see UaNum in
    // solidRegionDiffusionNo.H), same damped-growth treatment.
    scalar maxDeltaTUa = maxUa/(UaNum + small);

    // 2026-09-29: report which criterion limits the time step (a run whose
    // step stayed near 1e-5 s gave no clue in the log).
    Info<< "dtLimits Co " << maxDeltaTFluid << " Di " << maxDeltaTSolid
        << " Ra " << maxDeltaTRad << " Ga " << maxDeltaTGas
        << " Ua " << maxDeltaTUa << " (ratios dtMax/dt)" << endl;

    const scalar dt0 = runTime.deltaT().value() *
        min
        (
            min
            (
                min
                (
                    min
                    (
                        min
                        (
                            min(maxDeltaTFluid, 1.0 + 0.1 * maxDeltaTFluid),
                            min(maxDeltaTSolid, 1.0 + 0.1 * maxDeltaTSolid)
                        ),
                        min(maxDeltaTRad, 1.0 + 0.1 * maxDeltaTRad)
                    ),
                    min(maxDeltaTGas, 1.0 + 0.1 * maxDeltaTGas)
                ),
                min(maxDeltaTUa, 1.0 + 0.1 * maxDeltaTUa)
            ),
            1.2
        );

    runTime.setDeltaT
    (
        min(dt0,maxDeltaT)
    );
}
// ************************************************************************* //
