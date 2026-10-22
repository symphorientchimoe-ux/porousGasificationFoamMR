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
Class
    Foam::boundedArrheniusReactionRate

Description
    Arrhenius reaction rate given by:

        k = A * T^beta * exp(-Ta/T)

    Identical to the stock OpenFOAM ArrheniusReactionRate, except T is
    floored to a small positive value before it reaches pow(T, beta) or
    exp(-Ta/T). Gas-phase reactive ODE solvers (seulex, RKF45, ...) trial-
    evaluate this at intermediate/extrapolated states that are not
    physically bounded the way a converged field value is -- an unguarded
    T<=0 makes -Ta/T diverge (overflow for T<0, division-by-zero for
    T==0), trapped as SIGFPE under FOAM_SIGFPE. Flooring recovers the
    correct physical limit anyway (rate -> 0 as T -> 0+) instead of
    crashing. See porousGasificationFoam_OF8/gasifier_maillageReel_v6_
    the development notes for the root-cause history.

    Registered as reaction type "irreversibleBoundedArrhenius" (via
    makeBoundedArrheniusReactions.C), selectable in chemistryProperties'
    reactions{} exactly like the stock "irreversibleArrhenius".

SourceFiles
    boundedArrheniusReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef boundedArrheniusReactionRate_H
#define boundedArrheniusReactionRate_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class boundedArrheniusReactionRate;

Ostream& operator<<(Ostream&, const boundedArrheniusReactionRate&);


/*---------------------------------------------------------------------------*\
                class boundedArrheniusReactionRate Declaration
\*---------------------------------------------------------------------------*/

class boundedArrheniusReactionRate
{
    // Private Data

        scalar A_;
        scalar beta_;
        scalar Ta_;


public:

    // Constructors

        //- Construct from components
        inline boundedArrheniusReactionRate
        (
            const scalar A,
            const scalar beta,
            const scalar Ta
        );

        //- Construct from dictionary
        inline boundedArrheniusReactionRate
        (
            const speciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "boundedArrhenius";
        }

        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        inline scalar ddT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- Third-body efficiencies (beta = 1-alpha)
        //  non-empty only for third-body reactions
        //  with enhanced molecularity (alpha != 1)
        inline const List<Tuple2<label, scalar>>& beta() const;

        //- Species concentration derivative of the pressure dependent term
        //  By default this value is 1 as it multiplies the third-body term
        inline void dcidc
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalarField& dcidc
        ) const;

        //- Temperature derivative of the pressure dependent term
        //  By default this value is 0 since ddT of molecularity is approx.0
        inline scalar dcidT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const boundedArrheniusReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "boundedArrheniusReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
