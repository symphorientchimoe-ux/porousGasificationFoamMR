/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

Class
    Foam::solidPowerArrheniusReactionRate

Description
    Arrhenius reaction rate for solids with a temperature exponent :
        k = A T^n exp(-Ta/T)   for T >= Tcrit, 0 below.
    Input : (A n Ta Tcrit ...), the reaction then reads its own data.
    v1.2 (16/11/2026) : needed for the water resorption H2O -> moisture of
    Caraccio et al., Chem. Eng. J. 531 (2026) 174297 (n = 2, Ta < 0); also
    covers n = 0 (solidArrhenius) and n = 1 (solidTemperatureArrhenius).

SourceFiles
    solidPowerArrheniusReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef solidPowerArrheniusReactionRate_H
#define solidPowerArrheniusReactionRate_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of friend functions and operators
class solidPowerArrheniusReactionRate;
Ostream& operator<<(Ostream&, const solidPowerArrheniusReactionRate&);

/*---------------------------------------------------------------------------*\
                       Class solidPowerArrheniusReactionRate Declaration
\*---------------------------------------------------------------------------*/

class solidPowerArrheniusReactionRate
{
    // Private data

        scalar A_;
        scalar n_;
        scalar Ta_;
        scalar Tcrit_;


public:

    // Constructors

        //- Construct from components
        inline solidPowerArrheniusReactionRate
        (
            const scalar A,
            const scalar n,
            const scalar Ta,
            const scalar Tcrit
        );

        //- Construct from Istream
        inline solidPowerArrheniusReactionRate
        (
            const speciesTable& species,
            Istream& is
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "SolidPowerArrhenius";
        }

        inline scalar operator()
        (
            const scalar T,
            const scalar p,
            const scalarField& c
        ) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const solidPowerArrheniusReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidPowerArrheniusReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
