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
See strainMixingLimited.H for the model description.
\*---------------------------------------------------------------------------*/

#include "strainMixingLimited.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::strainMixingLimited<ReactionThermo>::strainMixingLimited
(
    const word& modelType,
    const ReactionThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    laminar<ReactionThermo>(modelType, thermo, turb, combustionProperties),
    Cmix_(this->coeffs().template lookupOrDefault<scalar>("Cmix", 1.0)),
    kappa_
    (
        IOobject
        (
            thermo.phasePropertyName(typeName + ":kappa"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    ),
    tmix_
    (
        IOobject
        (
            thermo.phasePropertyName(typeName + ":tmix"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimTime, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::strainMixingLimited<ReactionThermo>::~strainMixingLimited()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::strainMixingLimited<ReactionThermo>::correct()
{
    laminar<ReactionThermo>::correct();

    // U is not exposed by combustionModel/turbulence() (the solver owns it
    // and passes it explicitly to turbulence->divDevTau(U) in UEqn.H) --
    // looked up from the mesh registry instead, the standard OpenFOAM way
    // to reach a field not held by this class directly.
    const volVectorField& U =
        this->mesh().template lookupObject<volVectorField>("U");

    const volSymmTensorField Sij(symm(fvc::grad(U)));
    const volScalarField strainRate
    (
        sqrt(2.0)*mag(Sij) + dimensionedScalar("smallRate", dimless/dimTime, 1e-6)
    );

    tmix_ = Cmix_/strainRate;

    tmp<volScalarField> ttc(this->tc());
    const scalarField& tc = ttc();

    const scalarField& tmix = tmix_;

    forAll(tc, i)
    {
        if (tmix[i] > small)
        {
            kappa_[i] = tc[i]/(tc[i] + tmix[i]);
        }
        else
        {
            kappa_[i] = 1.0;
        }
    }

    // Debug temporaire (2026-09-23) : diagnostic de la raideur numerique
    // observee (tests 35/36, tres lent). A retirer une fois tranche.
    Info<< "DEBUG-strainMixingLimited"
        << " strainRate min/max = " << gMin(strainRate) << ", " << gMax(strainRate)
        << " tc min/max = " << gMin(tc) << ", " << gMax(tc)
        << " tmix min/max = " << gMin(tmix_) << ", " << gMax(tmix_)
        << " kappa min/max = " << gMin(kappa_) << ", " << gMax(kappa_)
        << endl;
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::strainMixingLimited<ReactionThermo>::R
(
    volScalarField& Y
) const
{
    return kappa_*laminar<ReactionThermo>::R(Y);
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::strainMixingLimited<ReactionThermo>::Qdot() const
{
    return volScalarField::New
    (
        this->thermo().phasePropertyName(typeName + ":Qdot"),
        kappa_*laminar<ReactionThermo>::Qdot()
    );
}


template<class ReactionThermo>
bool Foam::combustionModels::strainMixingLimited<ReactionThermo>::read()
{
    if (laminar<ReactionThermo>::read())
    {
        Cmix_ = this->coeffs().template lookupOrDefault<scalar>("Cmix", 1.0);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
