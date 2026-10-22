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
    Implementation of boundedArrheniusReactionRate -- see
    boundedArrheniusReactionRate.H for why T is floored here.
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::boundedArrheniusReactionRate::boundedArrheniusReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta)
{}


inline Foam::boundedArrheniusReactionRate::boundedArrheniusReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(dict.lookup<scalar>("A")),
    beta_(dict.lookup<scalar>("beta")),
    Ta_(dict.lookup<scalar>("Ta"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::boundedArrheniusReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    const scalar Tc = max(T, small);

    scalar ak = A_;

    if (mag(beta_) > vSmall)
    {
        ak *= pow(Tc, beta_);
    }

    if (mag(Ta_) > vSmall)
    {
        ak *= exp(-Ta_/Tc);
    }

    return ak;
}


inline Foam::scalar Foam::boundedArrheniusReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    const scalar Tc = max(T, small);

    scalar ak = A_;

    if (mag(beta_) > vSmall)
    {
        ak *= pow(Tc, beta_);
    }

    if (mag(Ta_) > vSmall)
    {
        ak *= exp(-Ta_/Tc);
    }

    return ak*(beta_+Ta_/Tc)/Tc;
}


inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::boundedArrheniusReactionRate::beta() const
{
    return NullObjectRef<List<Tuple2<label, scalar>>>();
}


inline void Foam::boundedArrheniusReactionRate::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcidc
) const
{}


inline Foam::scalar Foam::boundedArrheniusReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
}


inline void Foam::boundedArrheniusReactionRate::write(Ostream& os) const
{
    writeEntry(os, "A", A_);
    writeEntry(os, "beta", beta_);
    writeEntry(os, "Ta", Ta_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const boundedArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
