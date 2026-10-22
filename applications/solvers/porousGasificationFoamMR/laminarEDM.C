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
See laminarEDM.H for the model description.
\*---------------------------------------------------------------------------*/

#include "laminarEDM.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::laminarEDM<ReactionThermo>::readCoeffs()
{
    const auto& species = this->thermo().composition().species();

    A_ = this->coeffs().template lookupOrDefault<scalar>("A", 4.0);
    B_ = this->coeffs().template lookupOrDefault<scalar>("B", 0.5);
    oxidanti_ = species[this->coeffs().template lookup<word>("oxidant")];

    const PtrList<entry> reactionEntries(this->coeffs().lookup("reactions"));

    if (reactionEntries.size() != this->chemistryPtr_->nReaction())
    {
        FatalErrorInFunction
            << "laminarEDMCoeffs lists " << reactionEntries.size()
            << " reactions but chemistryProperties has "
            << this->chemistryPtr_->nReaction() << exit(FatalError);
    }

    fueli_.setSize(reactionEntries.size());
    s_.setSize(reactionEntries.size());
    productsi_.setSize(reactionEntries.size());

    forAll(reactionEntries, r)
    {
        const dictionary& dict = reactionEntries[r].dict();

        fueli_[r] = species[dict.lookup<word>("fuel")];
        s_[r] = dict.lookup<scalar>("s");

        const wordList products(dict.lookup("products"));
        productsi_[r].setSize(products.size());
        forAll(products, pi)
        {
            productsi_[r][pi] = species[products[pi]];
        }

        Info<< "    laminarEDM reaction " << r << " ("
            << reactionEntries[r].keyword() << "): fuel "
            << species[fueli_[r]] << ", s = " << s_[r]
            << ", products " << products << endl;
    }

    Info<< "    laminarEDM A = " << A_ << ", B = " << B_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::laminarEDM<ReactionThermo>::laminarEDM
(
    const word& modelType,
    const ReactionThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    laminar<ReactionThermo>(modelType, thermo, turb, combustionProperties),
    A_(4.0),
    B_(0.5),
    oxidanti_(-1),
    fueli_(),
    s_(),
    productsi_(),
    RR_(thermo.composition().Y().size()),
    kappa_()
{
    readCoeffs();

    forAll(RR_, i)
    {
        RR_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    "laminarEDM:RR." + thermo.composition().species()[i],
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );
    }

    kappa_.setSize(fueli_.size());
    forAll(kappa_, r)
    {
        kappa_.set
        (
            r,
            new volScalarField
            (
                IOobject
                (
                    "laminarEDM:kappa_"
                  + thermo.composition().species()[fueli_[r]],
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimless, 1)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::laminarEDM<ReactionThermo>::~laminarEDM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::laminarEDM<ReactionThermo>::correct()
{
    // laminar::correct() (gas ODE integration) deliberately not called --
    // see laminarEDM.H.

    const fvMesh& mesh = this->mesh();
    const scalar deltaT = mesh.time().deltaTValue();
    const PtrList<volScalarField>& Y = this->thermo().composition().Y();
    const label nR = fueli_.size();

    const volVectorField& U = mesh.template lookupObject<volVectorField>("U");
    const volScalarField strainRate(sqrt(2.0)*mag(symm(fvc::grad(U))));

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    // Realizability cap weighted by the gas volume fraction (2026-09-26, code
    // analysis): the gas stores porosityF*rho*Y per unit volume, so the amount
    // that can be consumed in one step is porosityF*rho*Y/deltaT, not rho*Y/deltaT.
    // Opt-in (laminarEDMCoeffs { porosityWeightedCap yes; }) to keep running cases unchanged.
    const bool epsCap =
        this->coeffs().template lookupOrDefault<Switch>("porosityWeightedCap", false);
    scalarField eps(rho.size(), 1.0);
    if (epsCap && mesh.foundObject<volScalarField>("porosityF"))
    {
        eps = mesh.lookupObject<volScalarField>("porosityF").primitiveField();
    }

    // Instantaneous kinetic fuel consumption rate per reaction (>= 0)
    PtrList<scalarField> Rkin(nR);
    for (label r = 0; r < nR; r++)
    {
        Rkin.set
        (
            r,
            new scalarField(-this->chemistryPtr_->calculateRR(r, fueli_[r]))
        );
    }

    label nReacting = 0;
    label nMixing = 0;
    label nCapped = 0;

    forAll(rho, celli)
    {
        const scalar YO2 = max(Y[oxidanti_][celli], 0.0);
        scalar O2consumed = 0;
        bool reacting = false;
        bool mixing = false;
        bool capped = false;

        for (label r = 0; r < nR; r++)
        {
            const scalar Rk = max(Rkin[r][celli], 0.0);

            if (Rk < vSmall)
            {
                kappa_[r][celli] = 1.0;
                continue;
            }
            reacting = true;

            const scalar Yf = max(Y[fueli_[r]][celli], 0.0);

            scalar YP = 0;
            forAll(productsi_[r], pi)
            {
                YP += max(Y[productsi_[r][pi]][celli], 0.0);
            }

            const scalar Redm =
                A_*rho[celli]*strainRate[celli]
               *min(min(Yf, YO2/s_[r]), B_*YP/(1.0 + s_[r]));

            const scalar Rreal = eps[celli]*rho[celli]*min(Yf, YO2/s_[r])/deltaT;

            scalar R = min(Rk, Redm);
            if (Redm < Rk) mixing = true;
            if (Rreal < R)
            {
                R = Rreal;
                capped = true;
            }

            kappa_[r][celli] = R/Rk;
            O2consumed += s_[r]*R;
        }

        // Reactions together must not consume more O2 than present
        const scalar O2available = eps[celli]*rho[celli]*YO2/deltaT;
        if (O2consumed > O2available && O2consumed > vSmall)
        {
            const scalar f = O2available/O2consumed;
            for (label r = 0; r < nR; r++)
            {
                kappa_[r][celli] *= f;
            }
            capped = true;
        }

        if (reacting) nReacting++;
        if (mixing) nMixing++;
        if (capped) nCapped++;
    }

    forAll(RR_, i)
    {
        RR_[i] = dimensionedScalar(dimMass/dimVolume/dimTime, 0);

        for (label r = 0; r < nR; r++)
        {
            RR_[i].field() +=
                kappa_[r].primitiveField()
               *this->chemistryPtr_->calculateRR(r, i)().field();
        }
    }

    forAll(kappa_, r)
    {
        kappa_[r].correctBoundaryConditions();
    }

    Info<< "laminarEDM: strainRate max = " << gMax(strainRate)
        << ", reacting cells = " << returnReduce(nReacting, sumOp<label>())
        << ", mixing-limited = " << returnReduce(nMixing, sumOp<label>())
        << ", realizability-capped = " << returnReduce(nCapped, sumOp<label>())
        << endl;
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::laminarEDM<ReactionThermo>::R
(
    volScalarField& Y
) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));
    fvScalarMatrix& Su = tSu.ref();

    const label specieI =
        this->thermo().composition().species()[Y.member()];
    Su += RR_[specieI];

    return tSu;
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::laminarEDM<ReactionThermo>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            this->thermo().phasePropertyName(typeName + ":Qdot"),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    scalarField& Qdot = tQdot.ref().primitiveFieldRef();

    forAll(RR_, i)
    {
        const scalar hi = this->thermo().composition().Hf(i);
        Qdot -= hi*RR_[i].field();
    }

    return tQdot;
}


template<class ReactionThermo>
bool Foam::combustionModels::laminarEDM<ReactionThermo>::read()
{
    if (laminar<ReactionThermo>::read())
    {
        A_ = this->coeffs().template lookupOrDefault<scalar>("A", 4.0);
        B_ = this->coeffs().template lookupOrDefault<scalar>("B", 0.5);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
