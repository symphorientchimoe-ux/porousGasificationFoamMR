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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
inline Foam::piecewiseLinearSolidTransport<thermo>::piecewiseLinearSolidTransport
(
    const thermo& t,
    const List<Tuple2<scalar, scalar>>& points
)
:
    thermo(t),
    points_(points)
{}


template<class thermo>
inline Foam::piecewiseLinearSolidTransport<thermo>::piecewiseLinearSolidTransport
(
    const word& name,
    const piecewiseLinearSolidTransport& pt
)
:
    thermo(name, pt),
    points_(pt.points_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::scalar Foam::piecewiseLinearSolidTransport<thermo>::K
(
    const scalar T
) const
{
    const label n = points_.size();

    if (n == 0)
    {
        return 0.0;
    }

    if (n == 1 || T <= points_[0].first())
    {
        return points_[0].second();
    }

    if (T >= points_[n-1].first())
    {
        return points_[n-1].second();
    }

    for (label i = 0; i < n-1; i++)
    {
        const scalar T0 = points_[i].first();
        const scalar T1 = points_[i+1].first();

        if (T >= T0 && T <= T1)
        {
            const scalar K0 = points_[i].second();
            const scalar K1 = points_[i+1].second();

            if (mag(T1 - T0) < SMALL)
            {
                return K0;
            }

            return K0 + (K1 - K0)*(T - T0)/(T1 - T0);
        }
    }

    // Should not reach here given the checks above
    return points_[n-1].second();
}


template<class thermo>
inline Foam::scalar Foam::piecewiseLinearSolidTransport<thermo>::alpha
(
    const scalar T
) const
{
    scalar Cp = this->Cp(T);

    scalar rho = this->rho(T);

    return (K(T)/(rho*Cp));
}


// ************************************************************************* //
