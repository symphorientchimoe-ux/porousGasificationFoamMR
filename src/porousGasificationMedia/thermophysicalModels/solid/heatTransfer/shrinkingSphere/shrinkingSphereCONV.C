/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 S.K. Tchimoe (porousGasificationFoamMR)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of porousGasificationFoamMR, derived from OpenFOAM.

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

#include "shrinkingSphereCONV.H"
#include "shrinkingParticle.H"
#include "Time.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(shrinkingSphereCONV, 0);
    addToRunTimeSelectionTable(heatTransferModel, shrinkingSphereCONV, porosity);
}

Foam::shrinkingSphereCONV::shrinkingSphereCONV
(
    const volScalarField& porosity,
    const volScalarField& initialPorosity
)
:
    heatTransferModel(porosity, initialPorosity),
    d0_(0),
    dMin_(0),
    correlated_(false),
    h_(0),
    Up_(db().lookupObject<volVectorField>("U")),
    rhop_(db().lookupObject<volScalarField>("rho")),
    alphap_(db().lookupObject<volScalarField>("thermo:alpha")),
    mup_(db().lookupObject<volScalarField>("thermo:mu")),
    thermop_(db().lookupObject<fluidThermo>("thermophysicalProperties"))
{
    read();
}

Foam::autoPtr<Foam::shrinkingSphereCONV> Foam::shrinkingSphereCONV::New
(
    const volScalarField& porosity,
    const volScalarField& initialPorosity
)
{
    return autoPtr<shrinkingSphereCONV>(new shrinkingSphereCONV(porosity, initialPorosity));
}

Foam::tmp<Foam::volScalarField> Foam::shrinkingSphereCONV::CONV() const
{
    tmp<volScalarField> tCONV
    (
        new volScalarField
        (
            IOobject
            (
                "CONVloc",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero", dimEnergy/dimTime/dimTemperature/dimVolume, 0.0
            )
        )
    );
    volScalarField& CONV = tCONV.ref();

    const scalarField& Cp = thermop_.Cp()();

    forAll(CONV, celli)
    {
        const scalar solid = max(1.0 - porosity()[celli], 0.0);
        const scalar phi0 =
            shrinkingParticle::phiS0(mesh_, initialPorosity(), celli);
        const scalar d =
            shrinkingParticle::diameter(d0_, dMin_, solid, phi0);
        const scalar av = shrinkingParticle::specificSurface(solid, d);

        scalar h = h_;
        if (correlated_)
        {
            const scalar Re = rhop_[celli]*mag(Up_[celli])*d/mup_[celli];
            const scalar kappa = Cp[celli]*alphap_[celli];
            const scalar Pr = mup_[celli]*Cp[celli]/max(kappa, small);
            h = (2.0 + 1.1*pow(Re, 0.6)*cbrt(Pr))*kappa/d;
        }
        CONV[celli] = h*av;
    }

    return tCONV;
}

bool Foam::shrinkingSphereCONV::read()
{
    IOdictionary dict
    (
        IOobject
        (
            "heatTransferProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    const dictionary& params = dict.subDict("Parameters");
    params.lookup("d0") >> d0_;
    dMin_ = params.lookupOrDefault<scalar>("dMin", 0.1*d0_);
    correlated_ = params.lookupOrDefault<word>("hModel", "constant") == "WakaoKaguei";
    if (!correlated_)
    {
        params.lookup("h") >> h_;
    }
    Info<< "shrinkingSphereCONV: d0 " << d0_ << " m, dMin " << dMin_
        << " m, h " << (correlated_ ? word("Wakao-Kaguei") : Foam::name(h_))
        << endl;
    return true;
}
