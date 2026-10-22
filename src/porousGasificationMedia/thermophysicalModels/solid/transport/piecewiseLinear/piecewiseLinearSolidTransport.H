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
    Foam::piecewiseLinearSolidTransport

Description
    Piecewise-linear K(T) properties for solid heat transport.
    Reads a list of (T K) breakpoints from the dictionary and linearly
    interpolates between the two bracketing points. Below the first
    point or above the last point, the boundary value is held constant
    (no extrapolation).

    Dictionary format:
        transport
        {
            breakpoints
            (
                (298    0.53)
                (1073   0.53)
                (1273   1.0)
            );
        }

    Templated into a given thermodynamics package.

SourceFiles
    piecewiseLinearSolidTransportI.H
    piecewiseLinearSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef piecewiseLinearSolidTransport_H
#define piecewiseLinearSolidTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Ostream;

// Forward declaration of friend functions and operators
template<class thermo> class piecewiseLinearSolidTransport;
template<class thermo>
Ostream& operator<<
(
    Ostream&,
    const piecewiseLinearSolidTransport<thermo>&
);


/*---------------------------------------------------------------------------*\
                  Class piecewiseLinearSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class thermo>
class piecewiseLinearSolidTransport
:
    public thermo
{
    // Private data

        //- (T, K) breakpoints, must be sorted by increasing T.
        List<Tuple2<scalar, scalar>> points_;

    // Private Member Functions

        //- Construct from components
        inline piecewiseLinearSolidTransport
        (
            const thermo& t,
            const List<Tuple2<scalar, scalar>>& points
        );


public:

    // Constructors

        //- Construct as named copy
        inline piecewiseLinearSolidTransport
        (
            const word&,
            const piecewiseLinearSolidTransport&
        );

        //- Construct from dictionary
        piecewiseLinearSolidTransport(const dictionary&);

    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "piecewiseLinear<" + thermo::typeName() + '>';
        }

        //- Thermal conductivity [W/mK]
        inline scalar K(const scalar T) const;

        //- Thermal diffusivity [kg/ms]
        inline scalar alpha(const scalar T) const;


    // Ostream Operator
        friend Ostream& operator<< <thermo>
        (
            Ostream&,
            const piecewiseLinearSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "piecewiseLinearSolidTransportI.H"

#ifdef NoRepository
#   include "piecewiseLinearSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
