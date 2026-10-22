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

\*---------------------------------------------------------------------------*/

#include "piecewiseLinearSolidTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
Foam::piecewiseLinearSolidTransport<thermo>::piecewiseLinearSolidTransport
(
    const dictionary& dict
)
:
    thermo(dict),
    points_()
{
    const dictionary& subDict = dict.subDict("transport");
    subDict.lookup("breakpoints") >> points_;

    if (points_.size() < 1)
    {
        FatalErrorInFunction
            << "piecewiseLinearSolidTransport: 'breakpoints' must contain "
            << "at least one (T K) pair."
            << exit(FatalError);
    }

    // Basic sanity check: breakpoints must be sorted by increasing T
    for (label i = 0; i < points_.size()-1; i++)
    {
        if (points_[i].first() >= points_[i+1].first())
        {
            FatalErrorInFunction
                << "piecewiseLinearSolidTransport: 'breakpoints' must be "
                << "sorted by strictly increasing temperature. "
                << "Found " << points_[i].first()
                << " followed by " << points_[i+1].first()
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os, const piecewiseLinearSolidTransport<thermo>& et
)
{
    operator<<(os, static_cast<const thermo&>(et));
    os << tab << et.points_;

    os.check
    (
        "Ostream& operator<<(Ostream& os, const piecewiseLinearSolidTransport& et)"
    );

    return os;
}


// ************************************************************************* //
