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

#include "shrinkingSphereST.H"
#include "shrinkingParticle.H"
#include "Time.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(shrinkingSphereST, 0);
    addToRunTimeSelectionTable(specieTransferModel, shrinkingSphereST, porosity);
}

Foam::shrinkingSphereST::shrinkingSphereST
(
    const volScalarField& por,
    const volScalarField& por0
)
:
    specieTransferModel(por, por0),
    d0_(0),
    dMin_(0),
    correlated_(false),
    h_(0),
    Dref_(0),
    Tref_(300)
{
    read();
}

Foam::autoPtr<Foam::shrinkingSphereST> Foam::shrinkingSphereST::New
(
    const volScalarField& por,
    const volScalarField& por0
)
{
    return autoPtr<shrinkingSphereST>(new shrinkingSphereST(por, por0));
}

Foam::tmp<Foam::volScalarField> Foam::shrinkingSphereST::ST() const
{
    tmp<volScalarField> tST
    (
        new volScalarField
        (
            IOobject
            (
                "STloc",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimTime, 0.0)
        )
    );
    volScalarField& ST = tST.ref();

    forAll(ST, celli)
    {
        const scalar solid = max(1.0 - por()[celli], 0.0);
        const scalar phi0 = shrinkingParticle::phiS0(mesh_, por0(), celli);
        const scalar d = shrinkingParticle::diameter(d0_, dMin_, solid, phi0);
        const scalar av = shrinkingParticle::specificSurface(solid, d);

        scalar hm = h_;
        if (correlated_)
        {
            const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
            const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
            const volScalarField& mu = mesh_.lookupObject<volScalarField>("thermo:mu");
            const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
            const scalar D = Dref_*pow(max(T[celli], 1.0)/Tref_, 1.75);
            const scalar Re = rho[celli]*mag(U[celli])*d/mu[celli];
            const scalar Sc = mu[celli]/(rho[celli]*D);
            hm = (2.0 + 1.1*pow(Re, 0.6)*cbrt(Sc))*D/d;
        }
        ST[celli] = hm*av;
    }

    return tST;
}

bool Foam::shrinkingSphereST::read()
{
    IOdictionary dict
    (
        IOobject
        (
            "specieTransferProperties",
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
    correlated_ = params.lookupOrDefault<word>("hModel", "constant") == "WakaoFunazkri";
    if (correlated_)
    {
        params.lookup("Dref") >> Dref_;
        Tref_ = params.lookupOrDefault<scalar>("Tref", 300);
    }
    else
    {
        params.lookup("h") >> h_;
    }
    Info<< "shrinkingSphereST: d0 " << d0_ << " m, dMin " << dMin_
        << " m, h_m " << (correlated_ ? word("Wakao-Funazkri") : Foam::name(h_))
        << endl;
    return true;
}
