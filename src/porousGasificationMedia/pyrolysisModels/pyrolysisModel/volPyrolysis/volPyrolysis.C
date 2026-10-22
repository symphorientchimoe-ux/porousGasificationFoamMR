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

#include "volFields.H"
#include "volPyrolysis.H"
#include "dimensionSet.H"
#include "addToRunTimeSelectionTable.H"
#include "mapDistribute.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceInterpolate.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcVolumeIntegrate.H"
#include "fvMatrices.H"
#include "fvCFD.H"
#include "DynamicList.H"
#include "IOmanip.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heterogeneousPyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(volPyrolysis, 0);

addToRunTimeSelectionTable(heterogeneousPyrolysisModel, volPyrolysis, noRadiation);
addToRunTimeSelectionTable(heterogeneousPyrolysisModel, volPyrolysis, radiation);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void volPyrolysis::readReactingOneDimControls()
{
    const dictionary& solution = this->solution().subDict("PIMPLE");
    solution.lookup("nNonOrthogonalCorrectors") >> nNonOrthCorr_;
    time_.controlDict().lookup("maxDi") >> maxDiff_;
    // 2026-09-17: optional, default 0.3; existing cases without the keyword
    // run unchanged.
    maxRad_ = time_.controlDict().lookupOrDefault<scalar>("maxRa", 0.3);
    // 2026-09-21 : same optional/default-0.3 pattern as maxRa above, see
    // gasRelaxNo().
    maxGas_ = time_.controlDict().lookupOrDefault<scalar>("maxGa", 0.3);
    equilibrium_.readIfPresent("equilibrium",coeffs_);
}

bool volPyrolysis::read()
{
    if (heterogeneousPyrolysisModel::read())
    {
        readReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}

bool volPyrolysis::read(const dictionary& dict)
{
    if (heterogeneousPyrolysisModel::read(dict))
    {
        readReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}

void volPyrolysis::solveSpeciesMass()
{
// eqZx2uHGn045
// eqZx2uHGn046

    if (debug)
    {
        Info<< "volPyrolysis:+solveSpeciesMass()" << endl;
    }

    if (active_)
    {
        volScalarField Mt(0.0 * Ym_[0]);
        volScalarField rhoLoc
        (
            max(rho_ * (1. - porosity_), dimensionedScalar("minRho",dimMass/dimVolume, SMALL))
        );

        for (label i = 0; i < Ys_.size(); ++i)
        {

            volScalarField& Yi = Ys_[i];
            Yi.ref() *= whereIs_;
            volScalarField sRhoSi = solidChemistry_->RRs(i);

            fvScalarMatrix YsEqn
            (
                fvm::ddt(rhoLoc,Yi)
             ==
                sRhoSi

            );

            YsEqn.relax();
            YsEqn.solve("Ys");

            Yi.max(0.0);
            Mt += Yi * rho_;

            Info<< "solid "<< Ys_[i].name()
                << " equation solved. Sources min/max   = " << gMin(sRhoSi)
                << ", " << gMax(sRhoSi)
                << "; values min Y = " << gMin(Ys_[i])
                << " max Y = " << gMax(Ys_[i]) << endl;
        }

        for (label i = 0; i < Ys_.size(); ++i)
        {
            Ym_[i] = whereIs_ * Ys_[i] * rho_;
            Ys_[i] = Ym_[i] / max(Mt,dimensionedScalar("minMass", dimMass/dimVolume, SMALL));
            Ym_[i] *= (1. - porosity_);
        }
    }
}


void volPyrolysis::solveEnergy()
{
// eqZx2uHGn047

    if (debug)
    {
        Info<< "volPyrolysis::solveEnergy()" << endl;
    }

    if (active_)
    {

        volTensorField composedK(K_ * (1 - porosity_) * anisotropyK_);
        radiationSh_ = radiation_;

        // Reintroduction fix (2026-09-14, see the development notes) :
        // radiationSh_ (= radiation_ = heterogeneousRadiationModel::
        // solidSh(), evaluated at the current/old T_) used to be added to
        // TEqn below as a fully explicit source. With heterogeneousP1,
        // its -sigma*Ts^4 boundary-loss part can reach ~1e8 W/m3 once
        // combustion is self-sustained -- applied explicitly, nothing
        // prevents it from removing more energy than a cell holds in one
        // timestep, driving T negative (root cause of the divergence
        // series this investigation started from). Fix: split
        // radiationSh_ back into its data part (independent of T_) and
        // its Ts^4 part (coefficient radiationRp_, from
        // heterogeneousRadiationModel::solidRp()), then apply the Ts^4
        // part implicitly via fvm::Sp -- same linearisation OpenFOAM's
        // own P1 model uses for the gas phase (see
        // heterogeneousRadiationModel::ST()). radiationRu below is
        // reconstructed algebraically from the already-computed
        // radiationSh_ + radiationRp_*T_^4, which is exactly equal to
        // solidSh()'s data part and avoids threading a second field
        // through the whole pyrolysis-model selection table.
        // Reintroduction fix, part 3 (2026-09-15, see
        // the development notes) : floor T_ at 0K before the T^4 here
        // (matching heterogeneousP1's own TsBounded, so this
        // reconstruction still exactly cancels the T^4 part down to the
        // T-independent G*coeff data term) and before the T^3 used below
        // as the implicit linearisation coefficient. Without this, once
        // T_ went slightly negative from an unrelated transient, the
        // *odd* power T^3 coefficient in fvm::Sp also went negative,
        // flipping the implicit radiative term from a stabilising sink
        // into a destabilising source and driving T further negative
        // (observed directly: T reached -1918K here vs -385K without
        // this fix) -- pow4()/pow3() cannot distinguish "very hot" from
        // "very negative", same class of bug as the gas-chemistry
        // exp(-Ta/T) issue already fixed with boundedArrheniusReactionRate.
        const volScalarField TBounded(max(T_, dimensionedScalar(dimTemperature, 0)));
        volScalarField radiationRu(radiationSh_ + radiationRp_ * pow4(TBounded));

        // T_ just before solve() -- used in the DIAG print below (2026-09-17,
        // see the development notes) so the printed trajectory shows the
        // value entering this solve, not just the post-solve result.
        const scalarField ToldDebug(T_.primitiveField());

        if (equilibrium_)
        {}
        else
        {
            volScalarField rhoCp
            (
                max
                (
                    rho_ * solidThermo_.Cp() * (1 - porosity_),
                    dimensionedScalar("minRhoCp",dimEnergy/dimTemperature/dimVolume,SMALL)
                )
            );

            // Simplistic immersed boundary for heat transport in solid phase.
            fvScalarMatrix TLap
            (
                fvm::laplacian(composedK, T_)
            );
 
            // Setting face fluxes on the border of porous media to 0.
            whereIs_.correctBoundaryConditions();
            surfaceScalarField  whereIsPatch  = fvc::interpolate(whereIs_);
            forAll(whereIsPatch,faceI)
            {
               if ( (whereIsPatch[faceI] > 0) and (whereIsPatch[faceI] < 1) )
               {
                   TLap.upper()[faceI] = 0.;
               }
            }
            forAll(whereIsPatch.boundaryField(),patchI)
            {
               if (isA<processorPolyPatch>(mesh_.boundaryMesh()[patchI]))
               {
                   forAll(whereIsPatch.boundaryField()[patchI],faceI)
                   {
                       if ( (whereIsPatch.boundaryField()[patchI][faceI] > 0) and (whereIsPatch.boundaryField()[patchI][faceI] < 1) )
                       {
                           TLap.boundaryCoeffs()[patchI][faceI] = 0;
                           TLap.internalCoeffs()[patchI][faceI] = 0;
                       }
                   }
               }
            }

            TLap.diag() = 0;
            TLap.negSumDiag();
            // Correct on orthogonal meshes
            // For non-orthogona meshes
            // TLap.source() = 0. cancels the non-orthogonal correction from the
            // divergence of composedK, which is spurious on the edges of porous
            // media and negiligble inside porous media in relevant cases.
            TLap.source() = 0.;

            // Reintroduction fix, part 2 (2026-09-15, see
            // the development notes) : same treatment as the radiative
            // sink above, applied to the solid<->gas convective coupling.
            // heatTransfer()() (= CONV_*(T_-Tgas)*whereIs_, CONV_ already
            // whereIs_-masked) used to be added here as a raw explicit
            // source -- the non-subintegrate branch of heatTransfer()
            // even carries a pre-existing code comment warning "This
            // works only for small CONV otherwise oscillations appear."
            // Confirmed exactly this symptom while validating the
            // radiative fix : a cell right on the mesh axis AND on the
            // bed's top surface layer (surfF_=1, so also subject to the
            // largest radiative border term) relaxed smoothly (no more
            // single-step blowup, radiative fix working) but to an
            // unphysical negative equilibrium temperature over ~1s
            // before recovering -- consistent with a still-explicit,
            // moderately-large (CONV~2300 W/m3/K here) convective term
            // allowing transient overshoot despite the radiative term
            // now being implicit. Same fvm::Sp linearisation as above :
            // CONV_*(T_-Tgas) split into an implicit CONV_*T_new sink and
            // an explicit CONV_*Tgas source.
            const volScalarField Tgas(gasThermo_.T());

            fvScalarMatrix TEqn
            (
                fvm::ddt(rhoCp, T_)
              - TLap
              + fvm::Sp(CONV_, T_)
            ==
                chemistrySh_ // eqZx2uHGn004, eqZx2uHGn017
              + CONV_ * Tgas // eqZx2uHGn005
              - heatUpGas_
              + radiationRu
              - fvm::Sp(radiationRp_ * pow3(TBounded), T_)
            );

            TEqn.relax();

            // Held igniter (2026-09-10): Dirichlet constraint in the igniter
            // cells until igniterDuration, then released. It must be applied
            // to the linear system (setValues() before solve()), not by
            // overwriting T_ after the solution: overwriting made the fluxes
            // of the equation inconsistent with the final value, injected a
            // spurious energy at every step and diverged. setValues() is what
            // fixedTemperatureConstraint::constrain() does internally; here the
            // equation is already in temperature, so no enthalpy conversion.
            if (mesh_.time().value() < igniterDuration_ and !igniterCells_.empty())
            {
                TEqn.setValues
                (
                    igniterCells_,
                    scalarField(igniterCells_.size(), igniterTemperature_)
                );
            }

            TEqn.solve();
        }

        scalar minTemp = GREAT;
        scalar maxTemp = -GREAT;
        scalar areThere = 0;
        label minCellI = -1;

        forAll(T_,cellI)
        {
            if (whereIs_[cellI] == 1 )
            {
                areThere = 1;
                if (T_[cellI] < minTemp)
                {
                    minTemp = T_[cellI];
                    minCellI = cellI;
                }
                if (T_[cellI] > maxTemp)
                {
                    maxTemp = T_[cellI];
                }
            }
        }
        scalar localMinTemp = minTemp;
        reduce(maxTemp, maxOp<scalar>());
        reduce(minTemp, minOp<scalar>());
        reduce(areThere, maxOp<scalar>());

        if (areThere == 1)
        {
            Info<< " pyrolysis min/max(T) = " << minTemp << ", " << maxTemp << endl;
        }
        else
        {
            Info<< " no solid phase " << endl;
        }

        // Diagnostic (2026-09-11): term-by-term breakdown of TEqn in the cell
        // of minimum temperature when it falls below 250 K, to identify the
        // term (chemistrySh_, heatTransfer, heatUpGas_, radiationSh_, rhoCp)
        // that drives the drift. Prints only; no effect on the solution.
        if (areThere == 1 and minTemp < 250 and minCellI >= 0
            and mag(localMinTemp - minTemp) < 1e-6)
        {
            scalar rhoCpLocal =
                max
                (
                    rho_[minCellI] * solidThermo_.Cp()()[minCellI]
                  * (1 - porosity_[minCellI]),
                    scalar(SMALL)
                );
            scalar TgasLocal = gasThermo_.T()[minCellI];
            scalar HTlocal = CONV_[minCellI] * (T_[minCellI] - TgasLocal);
            Pout<< "DIAG cell " << minCellI
                << " Told=" << ToldDebug[minCellI]
                << " T=" << T_[minCellI]
                << " Tgas=" << TgasLocal
                << " porF=" << porosity_[minCellI]
                << " porMax=" << porosityMax_[minCellI]
                << " whereIs=" << whereIs_[minCellI]
                << " rhoCp=" << rhoCpLocal
                << " CONV=" << CONV_[minCellI]
                << " chemistrySh=" << chemistrySh_[minCellI]
                << " heatTransfer(CONV*(T-Tgas))=" << HTlocal
                << " heatUpGas=" << heatUpGas_[minCellI]
                << " radiationSh=" << radiationSh_[minCellI]
                // 2026-09-16 (see the development notes) : radiationRu
                // and radiationRp_ already computed above (line ~195) as
                // exact inputs to TEqn's radiative terms this solve used --
                // printing the RAW data instead of reconstructing G after
                // the fact from a 5s checkpoint. radiationRu == coeff*G
                // (algebraic identity, see heterogeneousP1::calculate()),
                // so G = radiationRu*4*sigma/radiationRp_ can be recovered
                // offline (sigma=5.67e-8) without hardcoding it here.
                << " radiationRu=" << radiationRu[minCellI]
                << " radiationRp=" << radiationRp_[minCellI]
                << " deltaT=" << time_.deltaTValue()
                << endl;
        }
    }
}

void volPyrolysis::calculateMassTransfer()
{
    if (infoOutput_)
    {
        totalGasMassFlux_ = fvc::domainIntegrate(solidChemistry_->RRg());
        totalHeatRR_ = fvc::domainIntegrate(chemistrySh_);

        addedGasMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRg()) * time_.deltaT();
        lostSolidMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRs()) * time_.deltaT();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

volPyrolysis::volPyrolysis
(
    const word& modelType,
    const fvMesh& mesh,
    HGSSolidThermo& solidThermo,
    psiReactionThermo& gasThermo,
    volScalarField& whereIs
)
:
    heterogeneousPyrolysisModel(modelType, mesh),
    porosity_(whereIs),
    porosityArch_
    (
        IOobject
        (
            "porosityF0",
            time_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    STmodel_(specieTransferModel::New(porosity_,porosityArch_)),
    ST_
    (
        IOobject
        (
            "STvol",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimless/dimTime, 0.0)
    ),
    gasThermo_(gasThermo),
    Ygas_(gasThermo.composition().Y()),
    solidThermo_(solidThermo),
    solidChemistry_
    (
        BasicSolidChemistryModel<HGSSolidThermo>::New
        (
            solidThermo_,
            Ygas_,
            gasThermo.thermoName()
        )
    ),
    kappa_(solidThermo_.kappa()),
    K_(solidThermo_.K()),
    rho_(solidThermo_.rho()),
    rho0_
    (
        IOobject
        (
            "rhos0",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimVolume, 0.0)
    ),
    Ys_(solidThermo_.composition().Y()),
    Ym_(Ys_.size()),
    Msolid_(Ys_.size()),
    Msolidtotal_(
        IOobject
        (
            "Msolidtotal",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    T_(solidThermo_.T()),
    equilibrium_(false),
    subintegrateSwitch_(false),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    maxRad_(0.3),
    maxGas_(0.3),
    porosity0_(whereIs),
    porosityMax_
    (
        IOobject
        (
            "porosityMax",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("one", dimless, 1.0)
    ),
    voidFraction_(whereIs),
    radiation_(whereIs),
    radiationRp_(whereIs),
    phiHsGas_
    (
        IOobject
        (
            "phiHsGass",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),
    chemistrySh_
    (
        IOobject
        (
            "solidChemistrySh",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    porositySource_
    (
        IOobject
        (
            "porosityS",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    ),
    whereIs_
    (
        IOobject
        (
            "whereIs",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(1.0)
    ),
    whereIsNot_
    (
        IOobject
        (
            "whereIsNot",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0.0)
    ),
    heatUpGas_
    (
        IOobject
        (
            "heatUpGas",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    HTmodel_(heatTransferModel::New(porosity_,porosityArch_)),
    CONV_
    (
        IOobject
        (
            "CONVvol",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimEnergy/dimTime/dimTemperature/dimVolume , 0.0)

    ),
    radiationSh_
    (
        IOobject
        (
            "radiationSh",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    ),
    anisotropyK_
    (
        IOobject
        (
            "anisotropyK",
            time_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor("one", dimless, tensor(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    ),
    surfF_(whereIs_),
    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(dimensionedScalar("zero", dimMass/dimTime, 0.0)),
    igniterCells_(),
    igniterTemperature_(0),
    igniterDuration_(0),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0)),
    timeChem_(1.0)
{
    mesh.setFluxRequired(T_.name());

    ST_ = STmodel_->ST()();
    CONV_ = CONV();
    rho0_.ref() = rho_.ref();

    // Porosity ceiling (2026-09-11): the ash is inert (no reaction consumes
    // it), so its mass per unit cell volume is constant; it is taken here from
    // the initial state, and the volume it always occupies (at its own density,
    // ashDensity) bounds porosity_. Previously evolvePorosity() clamped at
    // 0.9999 -> 1.0 regardless of the ash and then wrote T_ = -1 K (a
    // sentinel) in cells that still held solid ash, which led to a slow drift
    // of the temperature towards negative values.
    {
        scalar ashDensity = coeffs().lookupOrDefault<scalar>("ashDensity", 650.0);

        label ashIdx = -1;
        forAll(Ys_, i)
        {
            if (Ys_[i].name() == "Yash")
            {
                ashIdx = i;
                break;
            }
        }

        if (ashIdx != -1)
        {
            // Masse de cendre par unite de volume TOTALE de cellule.
            volScalarField ashMassPerCV
            (
                Ys_[ashIdx] * rho_ * (1.0 - porosity_)
            );

            porosityMax_ = 1.0 - (ashMassPerCV / dimensionedScalar("ashRho", dimMass/dimVolume, ashDensity));
            porosityMax_.max(0.0);
            porosityMax_.min(1.0);

            Info<< "volPyrolysis: porosityMax (plafond du a la cendre "
                << "residuelle) min/max = " << gMin(porosityMax_)
                << ", " << gMax(porosityMax_) << endl;
        }
        else
        {
            WarningInFunction
                << "Espece 'Yash' introuvable dans Ys_ -- ceiling de "
                << "porosite desactive (porosityMax_=1 partout, "
                << "comportement identique a l'ancien code)." << endl;
        }
    }

    subintegrateSwitch_ = coeffs().lookupOrDefault("subintegrateHeatTransfer",false);

    // Held igniter (2026-09-10), see volPyrolysis.H.
    {
        word igniterCellSetName
        (
            coeffs().lookupOrDefault<word>("igniterCellSet", word::null)
        );

        igniterTemperature_ =
            coeffs().lookupOrDefault<scalar>("igniterTemperature", 0);
        igniterDuration_ =
            coeffs().lookupOrDefault<scalar>("igniterDuration", 0);

        if (igniterCellSetName != word::null)
        {
            cellSet igniterSet(mesh_, igniterCellSetName);
            igniterCells_ = igniterSet.toc();

            Info<< "volPyrolysis: igniter '" << igniterCellSetName
                << "' (" << returnReduce(igniterCells_.size(), sumOp<label>())
                << " cells total) held at " << igniterTemperature_
                << "K for " << igniterDuration_ << "s" << endl;
        }
    }

    forAll(Ys_, fieldI)
    {
        Ym_.set
        (
            fieldI,
            new volScalarField
            (
                IOobject
                (
                    Ys_[fieldI].name() + "m",
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimMass/dimVolume, 0.0)
            )
        );
        Ym_[fieldI].ref() = Ys_[fieldI].ref() * rho_.ref();

        // Msolid is used for renormalization of mass fractions of solid species.
        Msolid_.set
        (
            fieldI,
            new volScalarField
            (
                IOobject
                (
                    Ys_[fieldI].name() + "m2",
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimMass, 0.0)
            )
        );
        forAll(rho_,cellI)
        {
            Msolid_[fieldI][cellI] =
                Ys_[fieldI][cellI]
               * rho_[cellI]
               * mesh_.V()[cellI]
               * (1 - porosity_[cellI]);
        }
    }

    forAll(rho_,cellI)
    {
        if (porosity_[cellI] == 1)
        {
             whereIsNot_[cellI] = 1;
             whereIs_[cellI] = 0;
             rho0_[cellI] = 3.14;
        }
    }
    if (active_)
    {
        read();
    }

    forAll(whereIs_,cellI)
    {
        if (whereIs_[cellI] == 1)
        {
            bool surfC = false;
            forAll(mesh_.cellCells()[cellI],cellJ)
            {
                    if (whereIs_[mesh_.cellCells()[cellI][cellJ]] == 0)
                    {
                        surfC = true;
                    }
            }
            if (surfC)
            {
                surfF_[cellI] = 1;
            }
        }
    }

    forAll(rho_,cellI)
    {
        Msolidtotal_[cellI] = 0.0;

        for (label i = 0; i < Ys_.size(); ++i)
        {
             Msolid_[i][cellI]=
                Ys_[i][cellI]
               *rho_[cellI]
               *mesh_.V()[cellI]
               *(1 - porosity_[cellI]);

             Msolidtotal_[cellI] += Msolid_[i][cellI];
        }

        for (label i = 0; i < Ys_.size(); ++i)
        {
            if(Msolidtotal_[cellI] > 0.0)
            {
               Ys_[i][cellI] = (Msolid_[i][cellI] / (Msolidtotal_[cellI]));
            }
        }
    }

    // Bed collapse (2026-09-26) : opt-in, no-op when bedCollapse off.
    initBedCollapse();
}

volPyrolysis::volPyrolysis
(
    const word& modelType,
    const fvMesh& mesh,
    HGSSolidThermo& solidThermo,
    psiReactionThermo& gasThermo,
    volScalarField& whereIs,
    volScalarField& radiation,
    volScalarField& radiationRp
)
:
    heterogeneousPyrolysisModel(modelType, mesh),
    porosity_(whereIs),
    porosityArch_
    (
        IOobject
        (
            "porosityF0",
            time_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    STmodel_(specieTransferModel::New(porosity_,porosityArch_)),
    ST_
    (
        IOobject
        (
            "STvol",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimless/dimTime, 0.0)
    ),
    gasThermo_(gasThermo),
    Ygas_(gasThermo.composition().Y()),
    solidThermo_(solidThermo),
    solidChemistry_
    (
        BasicSolidChemistryModel<HGSSolidThermo>::New
        (
            solidThermo_,
            Ygas_,
            gasThermo.thermoName()
        )
    ),
    kappa_(solidThermo_.kappa()),
    K_(solidThermo_.K()),
    rho_(solidThermo_.rho()),
    rho0_
    (
        IOobject
        (
            "rhos0",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimVolume, 0.0)
    ),
    Ys_(solidThermo_.composition().Y()),
    Ym_(Ys_.size()),
    Msolid_(Ys_.size()),
    Msolidtotal_(
        IOobject
        (
            "Msolidtotal",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    T_(solidThermo_.T()),
    equilibrium_(false),
    subintegrateSwitch_(false),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    maxRad_(0.3),
    maxGas_(0.3),
    porosity0_(whereIs),
    porosityMax_
    (
        IOobject
        (
            "porosityMax",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("one", dimless, 1.0)
    ),
    voidFraction_(whereIs),
    radiation_(radiation),
    radiationRp_(radiationRp),
    phiHsGas_
    (
        IOobject
        (
            "phiHsGass",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),
    chemistrySh_
    (
        IOobject
        (
            "solidChemistrySh",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    porositySource_
    (
        IOobject
        (
            "porosityS",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    ),
    whereIs_
    (
        IOobject
        (
            "whereIs",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0.0)
    ),
    whereIsNot_
    (
        IOobject
        (
            "whereIsNot",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(1.0)
    ),
    heatUpGas_
    (
        IOobject
        (
            "heatUpGas",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    HTmodel_(heatTransferModel::New(porosity_,porosityArch_)),
    CONV_
    (
        IOobject
        (
            "CONVvol",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimEnergy/dimTime/dimTemperature/dimVolume , 0.0)
    ),
    radiationSh_
    (
        IOobject
        (
            "radiationSh",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    ),
    anisotropyK_
    (
        IOobject
        (
            "anisotropyK",
            time_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor("one", dimless, tensor(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    ),
    surfF_
    (
        IOobject
        (
            "surfF",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0.0)
    ),
    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(dimensionedScalar("zero", dimMass/dimTime, 0.0)),
    igniterCells_(),
    igniterTemperature_(0),
    igniterDuration_(0),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0)),
    timeChem_(1.0)
{
    mesh.setFluxRequired(T_.name());

    ST_ = STmodel_->ST()();
    CONV_ = CONV();
    rho0_.ref() = rho_.ref();

    // Porosity ceiling (2026-09-11): the ash is inert (no reaction consumes
    // it), so its mass per unit cell volume is constant; it is taken here from
    // the initial state, and the volume it always occupies (at its own density,
    // ashDensity) bounds porosity_. Previously evolvePorosity() clamped at
    // 0.9999 -> 1.0 regardless of the ash and then wrote T_ = -1 K (a
    // sentinel) in cells that still held solid ash, which led to a slow drift
    // of the temperature towards negative values.
    {
        scalar ashDensity = coeffs().lookupOrDefault<scalar>("ashDensity", 650.0);

        label ashIdx = -1;
        forAll(Ys_, i)
        {
            if (Ys_[i].name() == "Yash")
            {
                ashIdx = i;
                break;
            }
        }

        if (ashIdx != -1)
        {
            // Masse de cendre par unite de volume TOTALE de cellule.
            volScalarField ashMassPerCV
            (
                Ys_[ashIdx] * rho_ * (1.0 - porosity_)
            );

            porosityMax_ = 1.0 - (ashMassPerCV / dimensionedScalar("ashRho", dimMass/dimVolume, ashDensity));
            porosityMax_.max(0.0);
            porosityMax_.min(1.0);

            Info<< "volPyrolysis: porosityMax (plafond du a la cendre "
                << "residuelle) min/max = " << gMin(porosityMax_)
                << ", " << gMax(porosityMax_) << endl;
        }
        else
        {
            WarningInFunction
                << "Espece 'Yash' introuvable dans Ys_ -- ceiling de "
                << "porosite desactive (porosityMax_=1 partout, "
                << "comportement identique a l'ancien code)." << endl;
        }
    }

    subintegrateSwitch_ = coeffs().lookupOrDefault("subintegrateHeatTransfer",false);

    // Held igniter (2026-09-10), see volPyrolysis.H.
    {
        word igniterCellSetName
        (
            coeffs().lookupOrDefault<word>("igniterCellSet", word::null)
        );

        igniterTemperature_ =
            coeffs().lookupOrDefault<scalar>("igniterTemperature", 0);
        igniterDuration_ =
            coeffs().lookupOrDefault<scalar>("igniterDuration", 0);

        if (igniterCellSetName != word::null)
        {
            cellSet igniterSet(mesh_, igniterCellSetName);
            igniterCells_ = igniterSet.toc();

            Info<< "volPyrolysis: igniter '" << igniterCellSetName
                << "' (" << returnReduce(igniterCells_.size(), sumOp<label>())
                << " cells total) held at " << igniterTemperature_
                << "K for " << igniterDuration_ << "s" << endl;
        }
    }

    forAll(Ys_, fieldI)
    {
          Ym_.set
          (
                fieldI,
                new volScalarField
                (
                    IOobject
                    (
                        Ys_[fieldI].name() + "m",
                        time_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("zero",dimMass/dimVolume,0.0)
                )
          );

          Ym_[fieldI].ref() = Ys_[fieldI].ref() * rho_;
    }

    forAll(rho_,cellI)
    {
        if (porosity_[cellI] < 1.)
        {
             whereIsNot_[cellI] = 0.;
             whereIs_[cellI] = 1.;
        }
        else
        {
             whereIsNot_[cellI] = 1.;
             whereIs_[cellI] = 0.;
         rho0_[cellI] = 3.14;
        }
    }

    if (active_)
    {
        read();
    }

    forAll(whereIs_, cellI)
    {
        if (whereIs_[cellI] == 1)
        {
            bool surfC = false;
            forAll(mesh_.cellCells()[cellI],cellJ)
            {
                if (whereIs_[mesh_.cellCells()[cellI][cellJ]] == 0)
                {
                    surfC = true;
                }
            }
            if (surfC) surfF_[cellI] = 1;
        }
    }

    // Bed collapse (2026-09-26) : opt-in, no-op when bedCollapse off.
    initBedCollapse();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

volPyrolysis::~volPyrolysis()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
volScalarField& volPyrolysis::rhoConst() const
{
    return rho_;
}

volScalarField& volPyrolysis::rho()
{
    return rho_;
}

const volScalarField& volPyrolysis::T() const
{
    return T_;
}

const tmp<volScalarField> volPyrolysis::Cp() const
{
    return solidThermo_.Cp();
}


const volScalarField& volPyrolysis::kappa() const
{
    return kappa_;
}

const volScalarField& volPyrolysis::K() const
{
    return K_;
}

const volScalarField& volPyrolysis::surf() const
{
    return surfF_;
}

scalar volPyrolysis::maxDiff() const
{
    return maxDiff_;
}

scalar volPyrolysis::solidRegionDiffNo() const
{
    scalar NumCprho = 0.0;
    scalar DiNum = 0.0;

    if (mesh_.nInternalFaces() > 0)
    {
        surfaceScalarField KrhoCpbyDelta
        (
            mesh_.surfaceInterpolation::deltaCoeffs()
          * fvc::interpolate(K_)
        );

        surfaceScalarField Cprho
        (
            fvc::interpolate(Cp()*rho_)
        );

        NumCprho = max(Cprho.ref()).value();
        reduce(NumCprho, maxOp<scalar>());
        if (NumCprho != 0.)
        {
            DiNum = gMax(KrhoCpbyDelta.internalField())*time_.deltaTValue()/NumCprho;
        }
        else
        {
            DiNum = SMALL;
        }
    }

    return DiNum;
}

scalar volPyrolysis::radiativeRelaxNo() const
{
    // 2026-09-17 (see the development notes) : same role as
    // solidRegionDiffNo() above but for the radiative source in TEqn --
    // radiationSh_/radiationRp_ are computed from the PREVIOUS radiation.H
    // call (this control runs before radiation.H/evolve() each timestep,
    // same lagged-estimate principle as Courant number using the current
    // velocity field). Root cause of the borderAs=180+ divergence
    // (runs 8-23) : a cell hot enough + a high enough coeff (as+borderAs)
    // cools faster than deltaT can resolve, so the implicit linearisation
    // (built around Ts_old) overshoots well past the true radiative
    // equilibrium before the next step corrects it -- not a bug in G,
    // verified directly (2026-09-17).
    //
    // First version (2026-09-17, run23) used only the emission
    // linearisation slope (radiationRp_*Ts_old^3) -- this VANISHES once
    // Ts_old is floored to 0 (already crashed), silently disabling the
    // constraint exactly when the absorption side (coeff*G, unaffected by
    // the Ts floor) can still swing Ts hard. Fixed here: use the actual
    // net radiative source magnitude (mag(radiationSh_), captures
    // whichever of absorption/emission dominates, never vanishes just
    // because Ts floored) against a fixed reference temperature scale
    // (100K) rather than Ts_old itself -- bounds the potential
    // temperature swing from radiation alone to maxRa*100K per step.
    scalar RaNum = 0.0;

    if (active_ and !equilibrium_)
    {
        volScalarField rhoCp
        (
            max
            (
                rho_ * solidThermo_.Cp() * (1 - porosity_),
                dimensionedScalar("minRhoCp",dimEnergy/dimTemperature/dimVolume,SMALL)
            )
        );

        const dimensionedScalar Tref("Tref", dimTemperature, 100.0);

        RaNum =
            gMax((mag(radiationSh_) / (rhoCp * Tref))())
          * time_.deltaTValue();
    }

    return RaNum;
}

scalar volPyrolysis::maxRadDiff() const
{
    return maxRad_;
}

scalar volPyrolysis::gasRelaxNo()
{
    // 2026-09-21 : same role as radiativeRelaxNo() above but for the
    // explicit gas-phase source terms in EEqn.H (heatTransfer*porosityF +
    // heatUpGas_) -- these had no equivalent dt criterion, even though
    // the gas thermal mass (porosity_*rho_gas*Cp_gas) is far smaller
    // than the solid's (rho_gas ~0.2-1.2 kg/m3 vs solid ~1000s kg/m3),
    // making gas T much more prone to the same overshoot-past-
    // equilibrium failure already fixed for Ts via maxRa. Investigated
    // 2026-09-21 (gasifier_maillageReel_v8_wbtGeometry, tests 3-7): T
    // gaz repeatedly overshoots the JANAF ceiling (~2000-2090K -> 3500K)
    // at t~0.3s regardless of gravity/air BC/h-SAV/gas chemistry/
    // radiation/char oxidation (all isolated and confirmed innocent) --
    // heatUpGas/heatTransfer observed locally at several MW/m3 with a
    // tiny gas thermal mass, matching this exact failure signature.
    // heatTransfer() is non-const (recomputes each call, not cached like
    // radiationSh_), hence this function is not const either.
    scalar GaNum = 0.0;

    if (active_ and !equilibrium_)
    {
        volScalarField rhoCpGas
        (
            max
            (
                gasThermo_.rho() * gasThermo_.Cp() * porosity_,
                dimensionedScalar
                (
                    "minRhoCpGas",
                    dimEnergy/dimTemperature/dimVolume,
                    SMALL
                )
            )
        );

        const dimensionedScalar Tref("Tref", dimTemperature, 100.0);

        volScalarField gasSh(mag(heatTransfer() * porosity_ + heatUpGas_));

        GaNum =
            gMax((gasSh / (rhoCpGas * Tref))())
          * time_.deltaTValue();
    }

    return GaNum;
}

scalar volPyrolysis::maxGasDiff() const
{
    return maxGas_;
}

scalar volPyrolysis::maxTime() const
{
return timeChem_;
}

Switch volPyrolysis::equilibrium() const
{
    return equilibrium_;
}

void volPyrolysis::preEvolveRegion() {
    // Added for completeness
    heterogeneousPyrolysisModel::preEvolveRegion();

    // Iterates over every cell and sets cells containing solid phase
    // as reacting cell.
    forAll(T_, cellI)
    {
        if ( active_ && whereIs_[cellI] != 0)
        {
            solidChemistry_->setCellReacting(cellI, true);
        }
        else
        {
            solidChemistry_->setCellReacting(cellI, false);
        }
    }

}

void volPyrolysis::evolveRegion()
{
    voidFraction_ = porosity_;

    if (equilibrium_)
    {}
    else
    {
        CONV_ = CONV();
    }

    ST_ = STmodel_->ST()();

    timeChem_ = solidChemistry_->solve
    (
        time_.value() - time_.deltaTValue(),
        time_.deltaTValue()
    );

    heatUpGas_ = heatUpGasCalc()();

    chemistrySh_ = solidChemistry_->Sh()(); // eqZx2uHGn004

    solveSpeciesMass();

    evolvePorosity();

    // Bed collapse (2026-09-26, opt-in, see collapseBed()) : every
    // collapseInterval seconds of simulated time. Never entered when
    // bedCollapse is off (default) -> identical to the previous code.
    if (bedCollapse_ and active_)
    {
        const scalar idx =
            std::floor(time_.value()/collapseInterval_ + 1e-9);

        if (idx > lastCollapseIndex_)
        {
            lastCollapseIndex_ = idx;
            collapseBed();
        }
    }

    solidThermo_.correct(); // eqZx2uHGn046

    if (equilibrium_)
    {}
    else
    {
        for (int nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
        {
            solveEnergy();
        }
    }

    calculateMassTransfer();

    info();
}

void volPyrolysis::solveEnergyOnly()
{
    // Radiation<->Ts recoupling sub-iteration (2026-09-16, see
    // the development notes) : called from the top-level solver after
    // re-running radiation->correct() with the Ts_ this loop just
    // produced, so radiation_/radiationRp_ (references to the top-level
    // radiationF/radiationRpF fields) reflect the CURRENT Ts_ instead of
    // being frozen at the value from the start of the timestep. Chemistry/
    // mass transfer/porosity deliberately NOT repeated here -- those must
    // only be integrated once per timestep (see evolveRegion() above);
    // only the energy equation depends on radiation and needs re-solving.
    if (active_ and !equilibrium_)
    {
        for (int nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
        {
            solveEnergy();
        }
    }
}

void volPyrolysis::evolvePorosity()
{
    if (active_)
    {
        porositySource_ = solidChemistry_->RRpor(T_)();

        volScalarField& por = porosity_;

        fvScalarMatrix porosityEqn
        (
            fvm::ddt(por)
         ==
            porositySource_
        );

        porosityEqn.solve("porosity");

        Info<< "prosity equation solved. Sources min/max   = " << gMin(porositySource_)
            << ", " << gMax(porositySource_);

        Info<< "; values min Y = " << gMin(por)
            <<" max Y = " << gMax(por) << endl;

        forAll(porosity_,cellI)
        {
            // Porosity ceiling (2026-09-11, see the constructor): porosityMax_,
            // not 1.0, accounts for the residual ash, which stays solid. No
            // T_ = -1 sentinel any more: the cell stays in solveEnergy()
            // (whereIs_ excludes only boundary fluxes and diagnostics, not
            // ddt(rhoCp, T)), and the residual ash keeps a real heat capacity,
            // so the equation remains well posed.
            if (porosity_[cellI] > porosityMax_[cellI])
            {
                porosity_[cellI] = porosityMax_[cellI];
            }
            if (porosity_[cellI] < 0.0001)
            {
                porosity_[cellI] = 0.0;
            }
        }

        // whereIs_/whereIsNot_/surfF_ update factored into
        // updateSolidMasks() (2026-09-26) so that collapseBed() can reuse
        // it. Same per-cell operations, same order -> identical results.
        updateSolidMasks();
    }
    else
    {}
}

void volPyrolysis::updateSolidMasks()
{
        forAll(porosity_,cellI)
        {
            if (porosity_[cellI] < 1.0)
            {
                whereIs_[cellI] = 1.0;
                whereIsNot_[cellI] = 0.0;
            }
            else
            {
                whereIs_[cellI] = 0.0;
                whereIsNot_[cellI] = 1.0;
            }
        }

        surfF_=0;
        whereIs_.correctBoundaryConditions();
        surfaceScalarField  whereIsPatch  = fvc::interpolate(whereIs_);
        forAll(whereIsPatch,faceI)
        {
            if ( (whereIsPatch[faceI] > 0) and (whereIsPatch[faceI] < 1) )
            {
                if (whereIs_[mesh_.owner()[faceI]] == 1)
                {
                    surfF_[mesh_.owner()[faceI]] = 1;
                }
                if (whereIs_[mesh_.neighbour()[faceI]] == 1)
                {
                    surfF_[mesh_.neighbour()[faceI]] = 1;
                }
            }
        }
        forAll(whereIsPatch.boundaryField(),patchI)
        {
            if (isA<processorPolyPatch>(mesh_.boundaryMesh()[patchI]))
            {
                forAll(whereIsPatch.boundaryField()[patchI],faceI)
                {
                    if ( (whereIsPatch.boundaryField()[patchI][faceI] > 0) and (whereIsPatch.boundaryField()[patchI][faceI] < 1) )
                    {
                        if (whereIs_[mesh_.owner()[faceI + mesh_.boundaryMesh()[patchI].start()]] == 1)
                        {
                            surfF_[mesh_.owner()[faceI + mesh_.boundaryMesh()[patchI].start()]] = 1;
                        }
                    }
                }
            }
        }
}

// * * * * * * * * * * * * Bed collapse (2026-09-26)  * * * * * * * * * * * //
//
// Opt-in bed shrinkage/collapse, Yang et al. 2005 (Fuel 84:2026, eq. 7).
// Keywords (pyrolysisCoeffs, all optional) :
//   bedCollapse       off;          // master switch (default off)
//   shrinkA           0.8;          // a in h = dz[phiS0 - a(phiS0-phis)]/phiS0Ref
//   collapseInterval  20;           // [s] simulated time between collapses
//   bedZoneName       biomass_bed;  // cellZone holding the bed
//   phiS0Ref          <auto>;       // default: V-weighted mean of
//                                   //  (1-porosityF0) over the bed cells
//   phiS0Init         porosityF0;   // phiS0 init when no phiS0 file:
//                                   //  porosityF0 (1-porosityF0) or
//                                   //  porosityF (1-porosityF)
//   collapseMinPhis   1e-3;         // top-cell residue merged downwards
//
// The bed zone must be a structured set of vertical columns (axis z,
// columns identified by r = sqrt(x^2+y^2) of the cell centres, layers by
// z of the cell centres). dr/dz are deduced from the mesh.

void volPyrolysis::initBedCollapse()
{
    bedCollapse_ = coeffs().lookupOrDefault<Switch>("bedCollapse", false);

    if (!bedCollapse_)
    {
        return;
    }

    shrinkA_ = coeffs().lookupOrDefault<scalar>("shrinkA", 0.8);
    collapseInterval_ =
        coeffs().lookupOrDefault<scalar>("collapseInterval", 20.0);
    bedZoneName_ =
        coeffs().lookupOrDefault<word>("bedZoneName", "biomass_bed");
    collapseMinPhis_ =
        coeffs().lookupOrDefault<scalar>("collapseMinPhis", 1e-3);
    const word phiS0Init
    (
        coeffs().lookupOrDefault<word>("phiS0Init", "porosityF0")
    );

    if (collapseInterval_ <= 0)
    {
        FatalErrorInFunction
            << "collapseInterval must be > 0" << exit(FatalError);
    }
    if (phiS0Init != "porosityF0" and phiS0Init != "porosityF")
    {
        FatalErrorInFunction
            << "phiS0Init must be porosityF0 or porosityF, not "
            << phiS0Init << exit(FatalError);
    }

    const label zoneID = mesh_.cellZones().findZoneID(bedZoneName_);
    if (zoneID == -1)
    {
        FatalErrorInFunction
            << "bedCollapse : cellZone " << bedZoneName_ << " not found"
            << exit(FatalError);
    }

    bedCells_ = mesh_.cellZones()[zoneID];

    lastCollapseIndex_ =
        std::floor(time_.value()/collapseInterval_ + 1e-9);

    // --- Geometry : cell centre r,z and vertical extent of each cell
    const vectorField& C = mesh_.C().primitiveField();
    const pointField& pts = mesh_.points();
    const labelListList& cPts = mesh_.cellPoints();

    scalarField locR(bedCells_.size());
    scalarField locZ(bedCells_.size());
    scalarField locZlo(bedCells_.size());
    scalarField locZhi(bedCells_.size());
    scalar minDz = GREAT;

    forAll(bedCells_, bi)
    {
        const label celli = bedCells_[bi];
        locR[bi] = Foam::sqrt(sqr(C[celli].x()) + sqr(C[celli].y()));
        locZ[bi] = C[celli].z();
        scalar zlo = GREAT, zhi = -GREAT;
        forAll(cPts[celli], j)
        {
            zlo = min(zlo, pts[cPts[celli][j]].z());
            zhi = max(zhi, pts[cPts[celli][j]].z());
        }
        locZlo[bi] = zlo;
        locZhi[bi] = zhi;
        minDz = min(minDz, zhi - zlo);
    }
    reduce(minDz, minOp<scalar>());
    const scalar tol = 1e-3*minDz;

    // Every processor gets all (r,z) and clusters them identically
    List<scalarField> allR(Pstream::nProcs());
    List<scalarField> allZ(Pstream::nProcs());
    allR[Pstream::myProcNo()] = locR;
    allZ[Pstream::myProcNo()] = locZ;
    Pstream::gatherList(allR);
    Pstream::scatterList(allR);
    Pstream::gatherList(allZ);
    Pstream::scatterList(allZ);

    auto cluster = [tol](const List<scalarField>& all)
    {
        DynamicList<scalar> v;
        forAll(all, proci)
        {
            forAll(all[proci], i)
            {
                v.append(all[proci][i]);
            }
        }
        scalarField sv(v);
        Foam::sort(sv);
        DynamicList<scalar> u;
        forAll(sv, i)
        {
            if (u.empty() or sv[i] - u.last() > tol)
            {
                u.append(sv[i]);
            }
        }
        return scalarField(u);
    };

    const scalarField uR(cluster(allR));
    const scalarField uZ(cluster(allZ));
    const label nCol = uR.size();
    const label nLay = uZ.size();
    nBedSlots_ = nCol*nLay;

    auto findIdx = [tol](const scalarField& u, const scalar x)
    {
        label best = -1;
        scalar bestD = GREAT;
        forAll(u, i)
        {
            const scalar d = mag(x - u[i]);
            if (d < bestD)
            {
                bestD = d;
                best = i;
            }
        }
        return best;
    };

    bedSlot_.setSize(bedCells_.size());

    List<scalar> gZlo(nBedSlots_, 0.0);
    List<scalar> gZhi(nBedSlots_, 0.0);
    List<scalar> gV(nBedSlots_, 0.0);
    List<scalar> gCount(nBedSlots_, 0.0);

    forAll(bedCells_, bi)
    {
        const label c = findIdx(uR, locR[bi]);
        const label k = findIdx(uZ, locZ[bi]);
        const label s = c*nLay + k;
        bedSlot_[bi] = s;
        gZlo[s] += locZlo[bi];
        gZhi[s] += locZhi[bi];
        gV[s] += mesh_.V()[bedCells_[bi]];
        gCount[s] += 1.0;
    }

    Pstream::listCombineGather(gZlo, plusEqOp<scalar>());
    Pstream::listCombineScatter(gZlo);
    Pstream::listCombineGather(gZhi, plusEqOp<scalar>());
    Pstream::listCombineScatter(gZhi);
    Pstream::listCombineGather(gV, plusEqOp<scalar>());
    Pstream::listCombineScatter(gV);
    Pstream::listCombineGather(gCount, plusEqOp<scalar>());
    Pstream::listCombineScatter(gCount);

    slotZlo_ = gZlo;
    slotZhi_ = gZhi;
    slotV_ = gV;

    colSlots_.setSize(nCol);
    label nGap = 0;
    forAll(colSlots_, c)
    {
        DynamicList<label> sl;
        for (label k = 0; k < nLay; ++k)
        {
            const label s = c*nLay + k;
            if (gCount[s] > 1.5)
            {
                FatalErrorInFunction
                    << "bedCollapse : " << gCount[s] << " cells share the "
                    << "same (column,layer) slot (c=" << c << ", k=" << k
                    << ") -- the bed zone is not a structured wedge of "
                    << "vertical columns" << exit(FatalError);
            }
            if (gCount[s] > 0.5)
            {
                if (sl.size() and mag(slotZlo_[s] - slotZhi_[sl.last()]) > tol)
                {
                    ++nGap;
                }
                sl.append(s);
            }
        }
        colSlots_[c] = sl;
    }

    // --- Tracer phiS0
    IOobject hdr
    (
        "phiS0",
        time_.timeName(),
        mesh_,
        IOobject::NO_READ
    );
    const bool found = hdr.typeHeaderOk<volScalarField>(true);

    if (found)
    {
        phiS0Ptr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "phiS0",
                    time_.timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }
    else
    {
        phiS0Ptr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "phiS0",
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
        volScalarField& phiS0 = phiS0Ptr_();
        const volScalarField& porInit =
            (phiS0Init == "porosityF0") ? porosityArch_ : porosity_;
        forAll(bedCells_, bi)
        {
            const label celli = bedCells_[bi];
            phiS0[celli] = max(1.0 - porInit[celli], 0.0);
        }
        phiS0.correctBoundaryConditions();
    }

    // --- Reference solid fraction : from porosityF0 (restart-stable)
    {
        scalar sV = 0, sPV = 0;
        forAll(bedCells_, bi)
        {
            const label celli = bedCells_[bi];
            const scalar p0 = 1.0 - porosityArch_[celli];
            if (p0 > SMALL)
            {
                sV += mesh_.V()[celli];
                sPV += p0*mesh_.V()[celli];
            }
        }
        reduce(sV, sumOp<scalar>());
        reduce(sPV, sumOp<scalar>());
        const scalar autoRef = (sV > 0) ? sPV/sV : 1.0;
        phiS0Ref_ = coeffs().lookupOrDefault<scalar>("phiS0Ref", autoRef);
    }

    Info<< "volPyrolysis: bedCollapse ON (Yang 2005 eq.7) : zone "
        << bedZoneName_ << ", " << nCol << " columns x " << nLay
        << " layers (" << returnReduce(bedCells_.size(), sumOp<label>())
        << " cells, dz min " << minDz << " m, " << nGap << " gaps)"
        << ", shrinkA " << shrinkA_
        << ", collapseInterval " << collapseInterval_ << " s"
        << ", phiS0Ref " << phiS0Ref_
        << ", phiS0 " << (found ? "read from file" : "initialised from 1-")
        << (found ? word::null : phiS0Init)
        << ", collapseMinPhis " << collapseMinPhis_ << endl;
}


void volPyrolysis::collapseBed()
{
    const label nS = Ys_.size();
    const label iP = nS;        // phiS0 volume [m3]
    const label iH = nS + 1;    // solid sensible enthalpy [J]
    const label iVs = nS + 2;   // solid volume [m3]
    const label nQ = nS + 3;

    volScalarField& phiS0 = phiS0Ptr_();
    basicSolidMixture& mix = solidThermo_.composition();
    const scalarField& V = mesh_.V();
    const volScalarField& Tgas = gasThermo_.T();

    // rho_ consistent with the current Ys_/T_ (rho_ = harmonic mixture
    // density = sum(Y_i/rho_i)^-1), so that M = rho_(1-porosity)V and the
    // solid volume sum(M_i/rho_i) = (1-porosity)V hold exactly before the
    // remap. correct() is idempotent (recomputed again after collapseBed).
    solidThermo_.correct();

    // --- Pack per-slot quantities (global lists, one contributor per slot)
    List<scalar> q(nQ*nBedSlots_, 0.0);
    forAll(bedCells_, bi)
    {
        const label celli = bedCells_[bi];
        const label s = bedSlot_[bi];
        const scalar phis = 1.0 - porosity_[celli];
        const scalar Mc = rho_[celli]*phis*V[celli];
        for (label i = 0; i < nS; ++i)
        {
            q[s*nQ + i] = Ys_[i][celli]*Mc;
        }
        q[s*nQ + iP] = phiS0[celli]*V[celli];
        q[s*nQ + iH] = Mc*mix.hs(T_[celli], celli);
        q[s*nQ + iVs] = phis*V[celli];
    }
    Pstream::listCombineGather(q, plusEqOp<scalar>());
    Pstream::listCombineScatter(q);

    scalarField before(nQ, 0.0);
    for (label s = 0; s < nBedSlots_; ++s)
    {
        for (label j = 0; j < nQ; ++j)
        {
            before[j] += q[s*nQ + j];
        }
    }

    // --- Remap along each column
    List<scalar> qn(nQ*nBedSlots_, 0.0);
    label nClamped = 0;
    label nMerged = 0;
    scalar topMinBefore = GREAT, topMaxBefore = -GREAT;
    scalar topMinAfter = GREAT, topMaxAfter = -GREAT;
    scalar topAxisBefore = 0, topAxisAfter = 0;

    forAll(colSlots_, c)
    {
        const labelList& sl = colSlots_[c];
        if (sl.empty())
        {
            continue;
        }

        const scalar zBase = slotZlo_[sl[0]];
        scalar zTop = zBase;
        scalar topBefore = zBase;

        forAll(sl, kk)
        {
            const label s = sl[kk];
            const scalar dz = slotZhi_[s] - slotZlo_[s];
            const scalar Vk = slotV_[s];
            const scalar p0 = q[s*nQ + iP]/Vk;
            const scalar phis = q[s*nQ + iVs]/Vk;

            bool hasContent = false;
            for (label j = 0; j < nQ; ++j)
            {
                if (q[s*nQ + j] != 0)
                {
                    hasContent = true;
                }
            }
            if (!hasContent)
            {
                continue;
            }
            if (phis > 0)
            {
                topBefore = slotZhi_[s];
            }

            scalar h = dz*(p0 - shrinkA_*(p0 - phis))/phiS0Ref_;
            h = max(h, 0.0);
            if (h > dz)
            {
                if (h > dz*(1 + 1e-9))
                {
                    ++nClamped;
                }
                h = dz;
            }

            const scalar lo = zTop;
            const scalar hi = zTop + h;
            zTop = hi;

            // Targets overlapped by [lo,hi] ; the last one takes the
            // remainder so that the fractions sum to exactly 1.
            DynamicList<label> tgt;
            DynamicList<scalar> frac;
            if (h > 1e-12*dz)
            {
                forAll(sl, jj)
                {
                    const label t = sl[jj];
                    const scalar ov =
                        min(hi, slotZhi_[t]) - max(lo, slotZlo_[t]);
                    if (ov > 0)
                    {
                        tgt.append(t);
                        frac.append(ov/h);
                    }
                }
            }
            if (tgt.empty())
            {
                // zero height (or outside every layer) : whole content to
                // the layer containing lo
                label t = sl.last();
                forAll(sl, jj)
                {
                    if (lo < slotZhi_[sl[jj]])
                    {
                        t = sl[jj];
                        break;
                    }
                }
                tgt.append(t);
                frac.append(1.0);
            }

            scalar fSum = 0;
            forAll(tgt, m)
            {
                const scalar f =
                    (m == tgt.size() - 1) ? (1.0 - fSum) : frac[m];
                fSum += frac[m];
                for (label j = 0; j < nQ; ++j)
                {
                    qn[tgt[m]*nQ + j] += f*q[s*nQ + j];
                }
            }
        }

        // Thin residue in the new top cell -> merged into the cell below
        for (label kk = sl.size() - 1; kk > 0; --kk)
        {
            const label s = sl[kk];
            scalar Ms = 0;
            for (label i = 0; i < nS; ++i)
            {
                Ms += qn[s*nQ + i];
            }
            if (Ms <= 0 and qn[s*nQ + iVs] <= 0)
            {
                continue;
            }
            if (qn[s*nQ + iVs]/slotV_[s] < collapseMinPhis_)
            {
                const label sb = sl[kk - 1];
                for (label j = 0; j < nQ; ++j)
                {
                    qn[sb*nQ + j] += qn[s*nQ + j];
                    qn[s*nQ + j] = 0;
                }
                ++nMerged;
            }
            break;
        }

        topMinBefore = min(topMinBefore, topBefore);
        topMaxBefore = max(topMaxBefore, topBefore);
        topMinAfter = min(topMinAfter, zTop);
        topMaxAfter = max(topMaxAfter, zTop);
        if (c == 0)
        {
            topAxisBefore = topBefore;
            topAxisAfter = zTop;
        }
    }

    // --- Rebuild the local cells
    label nEmptied = 0;
    label nNegPor = 0;
    scalar maxVsErr = 0;
    scalarField porOld(porosity_.primitiveField());
    scalarField TPre(T_.primitiveField());

    forAll(bedCells_, bi)
    {
        const label celli = bedCells_[bi];
        const label s = bedSlot_[bi];

        scalar Mc = 0;
        for (label i = 0; i < nS; ++i)
        {
            Mc += qn[s*nQ + i];
        }

        phiS0[celli] = qn[s*nQ + iP]/V[celli];

        if (Mc <= 0)
        {
            // Above the new bed top : gas cell (as in the constructor)
            if (porosity_[celli] < 1.0)
            {
                ++nEmptied;
            }
            porosity_[celli] = 1.0;
            for (label i = 0; i < nS; ++i)
            {
                Ys_[i][celli] = 0.0;
            }
            T_[celli] = Tgas[celli];
            rho0_[celli] = 3.14;
        }
        else
        {
            for (label i = 0; i < nS; ++i)
            {
                Ys_[i][celli] = qn[s*nQ + i]/Mc;
            }

            // T from the remapped sensible enthalpy (Newton on the solid
            // mixture hs/Cp, which read the updated Ys_ of this cell)
            const scalar Hc = qn[s*nQ + iH];
            scalar T = T_[celli];
            for (label it = 0; it < 50; ++it)
            {
                const scalar f = Mc*mix.hs(T, celli) - Hc;
                const scalar dT = f/max(Mc*mix.Cp(T, celli), VSMALL);
                T -= dT;
                if (mag(dT) < 1e-10*max(T, 1.0))
                {
                    break;
                }
            }
            T_[celli] = T;

            // porosity from the component densities : mix.rho is the
            // harmonic mean, so Mc/rhoMix = sum_i M_i/rho_i(T)
            const scalar rhoMix = mix.rho(T, celli);
            scalar por = 1.0 - Mc/(rhoMix*V[celli]);
            maxVsErr = max
            (
                maxVsErr,
                mag((1.0 - por)*V[celli] - qn[s*nQ + iVs])/V[celli]
            );
            if (por < 0)
            {
                ++nNegPor;
                por = 0;
            }
            porosity_[celli] = por;
        }
    }

    // Gas side : keep the gas density unchanged through the jump of
    // porosity caused by the remap (fvm::ddt(porosityF, rho) in rhoEqn
    // would otherwise see it as a huge volume source/sink in one step).
    // Solid side : the remapped T is also the old-time value of this step
    // (fvm::ddt(rhoCp, T_) relaxes from T_.oldTime()).
    {
        scalarField& porO = porosity_.oldTime().primitiveFieldRef();
        scalarField& TO = T_.oldTime().primitiveFieldRef();
        const volScalarField rhoGas(gasThermo_.rho());
        scalar dGas = 0;
        forAll(bedCells_, bi)
        {
            const label celli = bedCells_[bi];
            const scalar dPor = porosity_[celli] - porOld[celli];
            porO[celli] = min(max(porO[celli] + dPor, 0.0), 1.0);
            TO[celli] = T_[celli];
            dGas += rhoGas[celli]*dPor*V[celli];
        }
        reduce(dGas, sumOp<scalar>());
        Info<< "bedCollapse: gas mass implicitly added by the porosity "
            << "jump (gas density kept) = " << dGas << " kg" << endl;
    }

    // Radiative source of this step (radiation.H, evaluated BEFORE the
    // collapse, at the pre-collapse Ts) : solveEnergy() rebuilds its
    // T-independent part as radiationSh + Rp*Ts^4 with the CURRENT Ts, so
    // a Ts changed by the remap would inject a spurious
    // Rp*(Ts_new^4 - Ts_old^4) (observed : emptied cells driven to
    // -1950 K in one step). Re-evaluate the emission at the remapped Ts
    // (same coefficients), and remove the radiative terms from the cells
    // that became gas. Only when radiation_ is a real radiation field
    // (the no-radiation constructor aliases it to the porosity field).
    if (&radiation_ != &porosity_ and &radiationRp_ != &porosity_)
    {
        forAll(bedCells_, bi)
        {
            const label celli = bedCells_[bi];
            if (porosity_[celli] >= 1.0)
            {
                radiation_[celli] = 0;
                radiationRp_[celli] = 0;
            }
            else
            {
                const scalar Tp = max(TPre[celli], 0.0);
                const scalar Tn = max(T_[celli], 0.0);
                radiation_[celli] +=
                    radiationRp_[celli]*(pow4(Tp) - pow4(Tn));
            }
        }
    }

    porosity_.correctBoundaryConditions();
    T_.correctBoundaryConditions();
    phiS0.correctBoundaryConditions();
    forAll(Ys_, i)
    {
        Ys_[i].correctBoundaryConditions();
    }

    // --- Everything that depends on porosity/Ys
    solidThermo_.correct();

    updateSolidMasks();

    forAll(Ys_, i)
    {
        Ym_[i] = whereIs_*Ys_[i]*rho_*(1.0 - porosity_);
    }

    // porosityMax_ : the ash moved (same formula as the constructor)
    {
        const scalar ashDensity =
            coeffs().lookupOrDefault<scalar>("ashDensity", 650.0);
        label ashIdx = -1;
        forAll(Ys_, i)
        {
            if (Ys_[i].name() == "Yash")
            {
                ashIdx = i;
            }
        }
        if (ashIdx != -1)
        {
            forAll(bedCells_, bi)
            {
                const label celli = bedCells_[bi];
                const scalar ashM =
                    Ys_[ashIdx][celli]*rho_[celli]*(1.0 - porosity_[celli]);
                porosityMax_[celli] =
                    min(max(1.0 - ashM/ashDensity, 0.0), 1.0);
            }
            porosityMax_.correctBoundaryConditions();
        }
    }

    if (!equilibrium_)
    {
        CONV_ = CONV();
    }
    heatUpGas_ = heatUpGasCalc()();
    chemistrySh_ *= whereIs_;

    // --- Conservation check (gSum over the cells, independent path)
    scalarField after(nQ, 0.0);
    forAll(bedCells_, bi)
    {
        const label celli = bedCells_[bi];
        const scalar Mc = rho_[celli]*(1.0 - porosity_[celli])*V[celli];
        for (label i = 0; i < nS; ++i)
        {
            after[i] += Ys_[i][celli]*Mc;
        }
        after[iP] += phiS0[celli]*V[celli];
        after[iH] += Mc*mix.hs(T_[celli], celli);
        after[iVs] += (1.0 - porosity_[celli])*V[celli];
    }
    forAll(after, j)
    {
        reduce(after[j], sumOp<scalar>());
    }

    Info<< "bedCollapse at t = " << time_.value()
        << " (shrinkA " << shrinkA_ << ", phiS0Ref " << phiS0Ref_ << ")"
        << nl
        << "    bed top axis column : " << topAxisBefore << " -> "
        << topAxisAfter << " m ; all columns min/max : "
        << topMinBefore << "/" << topMaxBefore << " -> "
        << topMinAfter << "/" << topMaxAfter << " m" << nl
        << "    cells emptied " << returnReduce(nEmptied, sumOp<label>())
        << ", top residues merged " << nMerged
        << ", h clamped to dz " << nClamped
        << ", porosity<0 clipped " << returnReduce(nNegPor, sumOp<label>())
        << ", max |Vs(rho_i) - Vs(remap)|/V "
        << returnReduce(maxVsErr, maxOp<scalar>()) << endl;

    for (label j = 0; j < nQ; ++j)
    {
        const string nm =
            (j < nS) ? string("M(" + Ys_[j].name() + ") [kg]")
          : (j == iP) ? string("phiS0 volume [m3]")
          : (j == iH) ? string("sensible enthalpy [J]")
          : string("solid volume [m3]");
        const scalar rel =
            (after[j] - before[j])/max(mag(before[j]), VSMALL);
        Info<< "    conservation " << nm.c_str() << " : before "
            << setprecision(15) << before[j] << " after " << after[j]
            << " rel " << setprecision(3) << rel << endl;
    }
    Info<< setprecision(IOstream::defaultPrecision());
}

Foam::tmp<Foam::volScalarField> volPyrolysis::Srho() const
{
    Foam::tmp<Foam::volScalarField> tSrho
    (
        new volScalarField
        (
            IOobject
            (
                "thermoSingleLayer::Srho",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    if (active_)
    {
        const speciesTable& gasTable = solidChemistry_->gasTable();

        forAll(gasTable,gasI)
        {
            tSrho = tSrho + solidChemistry_->RRg(gasI);;
        }
    }

    return tSrho;
}


Foam::tmp<Foam::volScalarField> volPyrolysis::Srho(const label i) const
{

    if (active_)
    {
        const speciesTable& gasTable = solidChemistry_->gasTable();

        label j = -1;
        forAll(gasTable,gasI)
        {
            if (gasTable[gasI] == Ygas_[i].name())
            {
                j = gasI;
            }
        }

        if (j > -1)
        {

            tmp<volScalarField> tRRiGas = solidChemistry_->RRg(j);   // 2026-09-29: upstream fix pjzuk c4ab0129 (gasTable index j, not Ygas_ index i)
            return tRRiGas;
        }
        else
        {
            return Foam::tmp<Foam::volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "volPyrolysis::Srho(" + Foam::name(i) + ")",
                        time_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_,
                    dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
                )
            );

        }
    }
    else
    {
        return Foam::tmp<Foam::volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "volPyrolysis::Srho(" + Foam::name(i) + ")",
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh_,
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            )
        );
    }
}

Foam::tmp<Foam::volScalarField> volPyrolysis::heatTransfer()
{
// eqZx2uHGn005
    Foam::tmp<Foam::volScalarField> Sh_ = Foam::tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "pyrolysisSh",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    if(equilibrium_)
    {}
    else
    {
        if (subintegrateSwitch_)
        {
            volScalarField rhoCpG(gasThermo_.rho() * gasThermo_.Cp() * porosity_);
            volScalarField Tgas = gasThermo_.T();
            volScalarField rhoCpS
                (
                        max
                        (
                            rho_ * solidThermo_.Cp() * (1 - porosity_),
                            dimensionedScalar("minRhoCp",dimEnergy/dimTemperature/dimVolume,SMALL)
                        )
                );
            volScalarField deltaTemp(T_ * 0);
            scalar deltaTime = time_.deltaTValue();

            forAll(deltaTemp,cellI)
            {
                if (whereIs_[cellI] == 1.0 && CONV_[cellI] > 0.0)
                {
                    deltaTemp[cellI] =
                        (Tgas[cellI]
                         * (rhoCpG[cellI]
                           +rhoCpS[cellI]
                           *exp(-(CONV_[cellI] * deltaTime*(rhoCpS[cellI]+rhoCpG[cellI]))
                                /(rhoCpS[cellI]*rhoCpG[cellI])))
                           +rhoCpS[cellI]*T_[cellI]
                           *(1 - exp( -(CONV_[cellI]*deltaTime*(rhoCpS[cellI]+rhoCpG[cellI]))/(rhoCpS[cellI]*rhoCpG[cellI])))
                        )
                        / (rhoCpS[cellI]+rhoCpG[cellI]) - Tgas[cellI];
                }
                else
                {
                    deltaTemp[cellI] = 0;
                }
            }

            forAll(Sh_(),cellI)
            {
                Sh_.ref()[cellI] = deltaTemp[cellI] * rhoCpG[cellI] * whereIs_[cellI] / deltaTime;
            }

            volScalarField HT(CONV_*(T_-Tgas));
            Info << "The heat transfer subintegration info:" << nl
                 << " no subintegration min, max:" << gMin(HT) << " " << gMax(HT) << nl
                 << " subintegration    min, max:" << gMin(Sh_()) << " " << gMax(Sh_()) << endl;
        }
        else
        {
            // This works only for small CONV otherwise oscillations appear.
            //
            // Reintroduction fix, part 5 (2026-09-15, see
            // the development notes) : floor T_ at 0K here too. This
            // field (Sh_, "heatTransfer") is consumed EXPLICITLY by the
            // GAS-phase energy equation (EEqn.H, "+ heatTransfer *
            // porosityF") in a different region than the one where T_ is
            // solved -- fixes 1/3/4 above made the SOLID-side TEqn robust
            // to a locally negative T_, but this cross-region explicit
            // term was untouched and still linear in raw T_, so once T_
            // dipped negative (from the radiative asymmetry at the
            // axis+top-surface corner cell, even bounded post-fix-4) this
            // term went to a magnitude of several MW/m3 and force-cooled
            // the local gas temperature Tgas, which then fed back into
            // the (now-implicit but Tgas-sourced) solid TEqn and stalled
            // its recovery at a deep unphysical equilibrium (observed:
            // T settling near -656K, mirroring a Tgas collapsed to
            // ~657K, instead of properly recovering the way run2's
            // ~-385K case slowly did with radiation-only fixes). Flooring
            // here prevents this specific field from ever carrying an
            // unphysical negative-T-amplified value into the gas region,
            // consistent with the same floor already applied on the
            // radiative side.
            const volScalarField TBounded(max(T_, dimensionedScalar(dimTemperature, 0)));
            volScalarField Tgas = gasThermo_.T();
            volScalarField HT(CONV_*(TBounded-Tgas));
            Sh_ = HT * whereIs_;
        }
    }
    return Sh_;
}

Foam::tmp<Foam::volScalarField> volPyrolysis::CONV() const
{
    Foam::tmp<Foam::volScalarField> CONVloc_ = Foam::tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "CONV",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimTemperature/dimTime/dimVolume, 0.0)
        )
    );

    if(equilibrium_)
    {}
    else
    {
        CONVloc_ = HTmodel_->CONV()*whereIs_;
    }
    return CONVloc_;
}

Foam::tmp<Foam::volScalarField> volPyrolysis::heatUpGasCalc() const
{

    Foam::tmp<Foam::volScalarField> hSh_ = Foam::tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "hSh_",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    if (active_)
    {
        volScalarField gasCp = gasThermo_.Cp();

        if (equilibrium_)
        {}
        else // eqZx2uHGn018
        {
            volScalarField tempSh = hSh_();
            // 2026-10-12 (energy balance, paper Section 3.7): heatUpGas switched off. The released gas leaves
            // the solid at Ts: the solid (equation in rho*Cp*dTs/dt) pays nothing more, and the gas receives
            // Srho(i)*hs_i(Ts) (term injSh of the solver). Previous form:
            // hSh_ = Cp_gas*(Ts - Tg)*Srho*whereIs*(1 - porosity), taken from the solid and given to the gas.
            tempSh = 0*gasThermo_.Cp() * (T_ - gasThermo_.T()) * Srho();
            hSh_ = tempSh * whereIs_ * (1 - porosity_);
        }
    }

    return hSh_;
}

Foam::tmp<Foam::volScalarField> volPyrolysis::heatUpGas() const
{
    return heatUpGas_;
}

Foam::tmp<Foam::volScalarField> volPyrolysis::solidChemistrySh() const
{
    return chemistrySh_;
}

void volPyrolysis::info() const
{
    Info<< "\nPyrolysis: " << endl;

    Info<< indent << "Total gas mass produced  [kg] = " << addedGasMass_.value() << nl
        << indent << "Total solid mass lost    [kg] = " << lostSolidMass_.value() << nl
        << indent << "Realese rate of pyrolysis gases  [kg/s] = " << totalGasMassFlux_.value() << nl
        << indent << "Heat release rate [J/s] = " << totalHeatRR_.value() << nl;

        if (timeChem_ < GREAT )
        {
            Info << indent << "Suggested chemical time step from heterogeneous reactions [s] = "
                 << timeChem_ << nl;
        }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace heterogeneousPyrolysisModels

// ************************************************************************* //
