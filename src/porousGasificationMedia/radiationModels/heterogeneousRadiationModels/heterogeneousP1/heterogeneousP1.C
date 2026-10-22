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

#include "heterogeneousP1.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "heterogeneousAbsorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"
#include "wedgePolyPatch.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
    defineTypeNameAndDebug(heterogeneousP1, 0);
    addToRunTimeSelectionTable
    (
        heterogeneousRadiationModel,
        heterogeneousP1,
        porosity
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousP1::heterogeneousP1
(
    const volScalarField& T,
    const volScalarField& porosityF,
    const volScalarField& surfF,
    const volScalarField& Ts
)
:
    // Fix (2026-09-17, see the development notes) : was
    // heterogeneousRadiationModel(typeName, T) -- the 2-arg base
    // constructor sets Ts_(T), silently binding the base class's solid-
    // temperature reference to the GAS temperature instead of the Ts
    // parameter received here. Root cause of the whole borderAs
    // investigation. Now uses the 3-arg constructor that correctly binds
    // Ts_(Ts).
    heterogeneousRadiationModel(typeName, T, Ts),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    qr_
    (
        // Written with the results and read back on restart (2026-10-22):
        // with qrRelax < 1 the wall flux is blended with its previous value,
        // so a restart that starts again from qr = 0 does not continue the
        // same trajectory. Cases without a qr file start from zero, as before.
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    as_
    (
        IOobject
        (
            "as",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    borderAs_
    (
        IOobject
        (
            "borderAs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    ),
    borderL_
    (
        dimensionedScalar("borderL", dimLength, 0.0)
    ),
    porosityF_(porosityF),
    surfFI_
    (
        surfF
    ),
    surfF_
    (
        IOobject
        (
            "surfF",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0.0)
    ),
    whereIs_
    (
        IOobject
        (
            "whereIs",
            mesh_.time().timeName(),
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
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0.0)
    ),
    solidSh_
    (
        IOobject
        (
            "radiationSolidSh",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    ),
    solidRp_
    (
        IOobject
        (
            "radiationSolidRp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime)/pow4(dimTemperature), 0)
    )
{
    forAll(porosityF_,cellI)
    {
        if (porosityF_[cellI] > (1.0 - pow(10.0, -8.0)))
        {
            whereIs_[cellI] = 0.0;
            whereIsNot_[cellI] = 1.0;
        }
        else
        {
            whereIs_[cellI] = 1.0;
            whereIsNot_[cellI] = 0.0;
        }

    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousP1::~heterogeneousP1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiationModels::heterogeneousP1::read()
{
    if (heterogeneousRadiationModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

void Foam::radiationModels::heterogeneousP1::calculate()
{
    a_ = heterogeneousAbsorptionEmission_->aCont();
    as_ = heterogeneousAbsorptionEmission_->asCont();
    borderAs_ = heterogeneousAbsorptionEmission_->borderAsCont();
    E_ = heterogeneousAbsorptionEmission_->ECont();
    borderL_  = heterogeneousAbsorptionEmission_->borderL();
    const volScalarField sigmaEff(scatter_->sigmaEff());

    volScalarField surfV = surfF_ * 0;

    if (((mesh_.geometricD()).x() + (mesh_.geometricD()).y() + (mesh_.geometricD()).z()) == 3)
    {
        surfV.ref() = borderL_ * pow(mesh_.V(), 2. / 3.) * surfFI_.internalField() * dimensionedScalar("tmp",dimensionSet(0, -3, 0, 0, 0),1.);
        surfF_.ref() = borderL_ / pow(mesh_.V(), 1. / 3.) * surfFI_.internalField();
    }
    else if (((mesh_.geometricD()).x() + (mesh_.geometricD()).y() + (mesh_.geometricD()).z()) == 1)
    {
        // Reintroduction fix (2026-09-14, see the development notes) :
        // the original formula here was
        // surfF_[cellI] = borderL/sqrt(cSurf/2), with cSurf the area of
        // the cell's wedge-empty-direction (front/back) faces. For an
        // axisymmetric wedge mesh those faces are triangular and their
        // area collapses to 0 as the cell approaches the symmetry axis
        // (x->0), making surfF_ -- and therefore the radiative loss
        // borderAs*surfF_*sigma*Ts^4 -- blow up there. Confirmed by
        // testing the implicit-coupling fix above: it correctly stopped
        // the catastrophic single-timestep runaway, but the cell that
        // then settled at an unphysical negative equilibrium temperature
        // (~-311K) was exactly an axis cell (x=0.0033m, the minimum x in
        // its processor). cSurf/the face-area geometry is a meshing
        // artifact of representing 3D axisymmetric geometry as a thin
        // wedge, not real physics -- so instead of patching cSurf with an
        // arbitrary floor, reuse the same non-singular, already-used-
        // elsewhere-in-this-model formula as the full-3D branch above
        // (borderL/V^(1/3), a generic cell length scale that is never
        // pathologically small the way a projected face area can be).
        //
        // Second fix (2026-09-25): V^(1/3)
        // is itself still wedge-biased. In a wedge V = theta*r*A_merid
        // (Pappus), so V^(1/3) ~ r^(1/3) and surfF_ ~ r^(-1/3): measured
        // on the bed-top row (z=0.374), surfF_ went 1.23 (axis cell 268,
        // x=0.0033) -> 0.36 (x=0.0975) for physically identical cells,
        // and cell 268 is exactly where G grew x20 before the t=420.16s
        // crash (test54). Use instead the meridional length
        // sqrt(V/(theta*r)) = sqrt(A_merid), independent of the wedge
        // angle and of r (0.0025-0.0029 m on that row, surfF_ 0.52-0.60).
        // Original V^(1/3) version kept in
        // heterogeneousP1.C.bak_2026-09-25_preMeridionalL.
        vector wedgeAxis(Zero);
        scalar wedgeAngle = 0;
        forAll(mesh_.boundaryMesh(), patchI)
        {
            if (isA<wedgePolyPatch>(mesh_.boundaryMesh()[patchI]))
            {
                const wedgePolyPatch& wpp =
                    refCast<const wedgePolyPatch>(mesh_.boundaryMesh()[patchI]);
                wedgeAxis = wpp.axis();
                wedgeAngle = 2.0*acos(min(wpp.cosAngle(), 1.0));
                break;
            }
        }

        if (wedgeAngle < small)
        {
            FatalErrorInFunction
                << "wedge geometry branch but no valid wedge patch found"
                << exit(FatalError);
        }

        forAll(surfFI_,cellI)
        {
            if(surfFI_[cellI] > 0)
            {
                const vector& c = mesh_.C()[cellI];
                const scalar r =
                    max(mag(c - (c & wedgeAxis)*wedgeAxis), small);
                surfF_[cellI] =
                    borderL_.value()
                   /sqrt(mesh_.V()[cellI]/(wedgeAngle*r));
            }
        }
    }
    else
    {
        vector f(0, 0 , 0);
        f.x() = max((mesh_.geometricD()).x(), 0);
        f.y() = max((mesh_.geometricD()).y(), 0);
        f.z() = max((mesh_.geometricD()).z(), 0);
        forAll(surfFI_,cellI)
        {
            if(surfFI_[cellI] > 0)
            {
                scalar cSurf = 0;
                forAll(mesh_.cells()[cellI],faceI)
                {
                    cSurf += mag(( f & (mesh_.Sf()[mesh_.cells()[cellI][faceI]]) ));
                }
                surfF_[cellI] = borderL_.value() / (2 * mesh_.V()[cellI] / cSurf);
            }
        }
    }

    scalar totalSurf = gSum(surfV);

    scalar totalVol = 0;
    forAll(porosityF_, cellI)
    {
        if (porosityF_[cellI] > (1.0 - pow(10.0, -8.0)))
        {
            whereIs_[cellI] = 0.0;
            whereIsNot_[cellI] = 1.0;
        }
        else
        {
            whereIs_[cellI] = 1.0;
            whereIsNot_[cellI] = 0.0;
            totalVol += mesh_.V()[cellI];
        }
    }

    reduce(totalVol, sumOp<scalar>());

    Info << "Radiation active volume to porous media volume ratio: "
         << totalSurf/max(totalVol,SMALL) << endl;

    // Construct diffusion
    const volScalarField gamma
    (
        IOobject
        (
            "gammaRad",
            G_.mesh().time().timeName(),
            G_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        1.0/(3.0*(a_ + as_*whereIs_ + borderAs_*surfF_ ) + sigmaEff)  // eqZx2uHGn010
    );

    // Reintroduction fix, part 3 (2026-09-15, see the development notes) :
    // floor Ts_ at 0K before raising to an even power here. sigma*Ts^4 is
    // a blackbody EMISSION term -- physically it must grow with how hot a
    // cell is, never with how far below absolute zero it (unphysically)
    // reads. But pow4() cannot tell "very hot" from "very negative" apart
    // (same sign, same huge magnitude) -- confirmed directly on run3
    // (radiationFixTest test case): once cell 418 (axis + bed-top-surface
    // cell) went slightly negative from an unrelated transient, this term
    // and its implicit linearisation coefficient in
    // volPyrolysis::solveEnergy() (fvm::Sp(solidRp*Ts^3, Ts), odd power,
    // sign-sensitive) flipped from a damping sink into a destabilising
    // one, driving Ts to -1918K (worse than the -385K seen without this
    // floor). Same principle as the pre-existing boundedArrheniusReactionRate
    // fix for gas chemistry (floor T before it enters an unguarded power/
    // exponential of T). solidRp_ below is recomputed from the SAME
    // floored temperature so volPyrolysis's algebraic reconstruction
    // (radiationRu = solidSh_ + solidRp_*T_old^4, see solveEnergy())
    // stays exactly consistent.
    const volScalarField TsBounded(max(Ts_, dimensionedScalar(dimTemperature, 0)));

    // eqZx2uHGn013
    volScalarField solidRadiation = (as_ * whereIs_ + borderAs_ * surfF_) * physicoChemical::sigma * pow4(TsBounded);

    // Solve G transport equation
    // eqZx2uHGn009
    solve
    (
        fvm::laplacian(gamma, G_)
        - fvm::Sp((a_ + as_ * whereIs_ + borderAs_ * surfF_ ), G_)
       ==
        - 4.0 * (a_ * physicoChemical::sigma * pow4(T_) + solidRadiation)
    );

    // DEBUG-radCoeff instrumentation (2026-09-25, test54) removed after
    // the fix was validated (test55) -- kept in
    // heterogeneousP1.C.bak_2026-09-25_withDebugRadCoeff.

    // Floor (2026-09-16, see the development notes) : G is a radiative
    // intensity integral (4*sigma*T^4 at equilibrium) and is physically
    // never negative. Kept as a cheap defensive floor even though the
    // 2026-09-17 investigation found the real cause of the borderAs
    // divergence elsewhere (Ts_ was bound to the gas temperature instead
    // of the solid one, see the 3-arg heterogeneousRadiationModel
    // constructor) -- G itself was never actually the problem once that
    // was fixed, but this floor costs nothing and the underlying
    // non-negativity is a genuine physical constraint regardless.
    G_ = max(G_, dimensionedScalar(G_.dimensions(), 0));

    // eqZx2uHGn012
    forAll(G_,cellI)
    {
        if (surfF_[cellI] == 0)
        {
            solidSh_[cellI] = (G_[cellI] * as_[cellI] - 4.0 * solidRadiation[cellI]) * whereIs_[cellI];
        }
        else
        {
            solidSh_[cellI] = (G_[cellI] * (as_[cellI] + surfF_[cellI] * borderAs_[cellI]) - 4.0 * solidRadiation[cellI]) * whereIs_[cellI];        }

        // Coefficient of Ts^4 in solidSh_ above (same combination as
        // solidRadiation, whereIs-masked the same way solidSh_ is) --
        // exposed via solidRp() so volPyrolysis::solveEnergy() can apply
        // this sink implicitly (fvm::Sp(solidRp*Ts^3, Ts)) instead of as
        // a raw explicit source. Reintroduction fix, 2026-09-14 -- see
        // the development notes.
        solidRp_[cellI] =
            4.0
          * (as_[cellI] * whereIs_[cellI] + borderAs_[cellI] * surfF_[cellI])
          * physicoChemical::sigma.value()
          * whereIs_[cellI];
    }

    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();

    // qr under-relaxation (2026-10-14) : the
    // coupled wall BC (turbulentTemperatureRadCoupledMixed) sets the face T
    // from qr, and qr comes from the Marshak BC on G, which uses that same
    // face T (4*sigma*Tp^4), one time step behind. This explicit loop has
    // gain ~ Ep*16*sigma*Tp^3/(kappa_s*delta_s + kappa_g*delta_f), close to
    // 1 on the freeboard wall; the bed collapse pushes it past 1 and the
    // face T alternates step to step with growing amplitude until < 0 K
    // (7 NaN in ign1000/M5). qrRelax (radiationProperties, default 1 =
    // unchanged behaviour) blends the new qr with the previous step's;
    // stable for qrRelax < 2/(1 + gain). On the first call (qr still 0)
    // the new value is taken as is.
    const scalar qrRelax = lookupOrDefault<scalar>("qrRelax", 1.0);

    // Calculate radiative heat flux on boundaries.
    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (!G_.boundaryField()[patchi].coupled())
        {
            const scalarField qrNew
            (
                -gamma.boundaryField()[patchi]
                *G_.boundaryField()[patchi].snGrad()
            );

            if (qrRelax < 1.0 && gMax(mag(qrBf[patchi])) > 0)
            {
                qrBf[patchi] =
                    qrRelax*qrNew + (1.0 - qrRelax)*qrBf[patchi];
            }
            else
            {
                qrBf[patchi] = qrNew;
            }
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::radiationModels::heterogeneousP1::Rp() const
{
    return volScalarField::New
    (
        "Rp",
        4.0*heterogeneousAbsorptionEmission().eCont()*physicoChemical::sigma
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiationModels::heterogeneousP1::Ru() const
{
    const volScalarField::Internal& G = G_();
    const volScalarField::Internal E = heterogeneousAbsorptionEmission_->ECont()()();
    const volScalarField::Internal a = heterogeneousAbsorptionEmission_->aCont()()();

    return a*G - E;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousP1::solidSh() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "tSolidSh",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            solidSh_
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousP1::solidRp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "tSolidRp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            solidRp_
        )
    );
}


// ************************************************************************* //
