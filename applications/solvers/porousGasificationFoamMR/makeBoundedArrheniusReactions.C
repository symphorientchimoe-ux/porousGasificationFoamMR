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

/*---------------------------------------------------------------------------*\
Description
    Registers IrreversibleReaction<ThermoPhysics, boundedArrheniusReactionRate>
    under runtime-selection type name "irreversibleBoundedArrhenius", for
    exactly the one ThermoPhysics combination this project's gas mixture
    uses (transport=sutherland, thermo=janaf, energy=sensibleEnthalpy,
    equationOfState=perfectGas, specie=specie -- matches
    thermophysicalProperties' thermoType and thermo.compressibleGas).

    Deliberately narrower than OpenFOAM's own makeReactions.C (which sweeps
    forCommonGases/forCommonLiquids/forPolynomials for every stock
    ReactionRate) -- only IrreversibleReaction is needed here (the only
    reaction type this project's gas mechanism uses), for only this one
    ThermoPhysics, to keep this file's scope minimal and its compile time
    negligible. See boundedArrheniusReactionRate.H for the fix itself.

    NOT a defineReaction -- Reaction<ThermoPhysics>'s runtime-selection
    table for this ThermoPhysics is already defined once by OpenFOAM's own
    makeReactions.C (linked via -lchemistryModel); redefining it here would
    be a duplicate-definition link error. Only the new ReactionRate variant
    is added to that existing table via addToRunTimeSelectionTable, exactly
    like OpenFOAM's own forCommonGases(makeIRNReactions, ArrheniusReactionRate)
    call does for the stock "irreversibleArrhenius" type.
\*---------------------------------------------------------------------------*/

#include "makeReaction.H"

#include "boundedArrheniusReactionRate.H"

#include "forCommonGases.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forThermo
    (
        sutherlandTransport,
        sensibleEnthalpy,
        janafThermo,
        perfectGas,
        specie,
        makeIReactions,
        boundedArrheniusReactionRate
    );
}

// ************************************************************************* //
