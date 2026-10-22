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
    Foam::heterogeneousPyrolysisModels::volPyrolysis

Description

    Reacting and radiating volumetric pyrolysis model
SourceFiles
    volPyrolysis.C

\*---------------------------------------------------------------------------*/

#ifndef volPyrolysis_H
#define volPyrolysis_H

#include "heterogeneousPyrolysisModel.H"
#include "HGSSolidThermo.H"
#include "BasicSolidChemistryModel.H"
#include "cellSet.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heterogeneousPyrolysisModels
{


/*---------------------------------------------------------------------------*\
                      Class volPyrolysis Declaration
\*---------------------------------------------------------------------------*/

class volPyrolysis
:
    public heterogeneousPyrolysisModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        volPyrolysis(const volPyrolysis&);

        //- Disallow default bitwise assignment
        void operator=(const volPyrolysis&);

        //- Read model controls
        void readReactingOneDimControls();


protected:

    //- Porosity
    volScalarField& porosity_;

    volScalarField porosityArch_;

    //- Reference to specie transport model
    autoPtr<specieTransferModel> STmodel_;

    volScalarField ST_;

    //- Thermophysical properties of gaseous species.
    psiReactionThermo& gasThermo_;

    //- List of gas components
    PtrList<volScalarField>& Ygas_;

    //- Thermophysical properties of solid phase
    HGSSolidThermo& solidThermo_;

    //- Solid chemistry model
    autoPtr<BasicSolidChemistryModel<HGSSolidThermo>> solidChemistry_;

    // Reference to solid thermo properties

        //- Absorption coefficient [1/m]
        const volScalarField& kappa_;

        //- Thermal conductivity [W/m/K]
        const volScalarField& K_;

        //- Density [kg/m3]
        volScalarField& rho_;

        //- Density for the initial moment
        volScalarField rho0_;

        //- List of solid components
        PtrList<volScalarField>& Ys_;

       //- List of solid component masses
        PtrList<volScalarField> Ym_;

        //- List of solid component masses [kg]
        PtrList<volScalarField> Msolid_;

        //- solid total mass [kg]
        volScalarField Msolidtotal_;

        // Non-const access to temperature
        volScalarField& T_;

    //- Name of the radiative flux in the primary region
    word primaryRadFluxName_;

    ///-Switch for equilibrium/nonequilibrium model
    Switch equilibrium_;

    //- Switch for temporal subintegration of heatTransfer model
    Switch subintegrateSwitch_;

    // Solution parameters

        //- Number of non-orthogonal correctors
        label nNonOrthCorr_;
        //- Maximum diffussivity
        scalar maxDiff_;
        //- Maximum radiative relaxation number (2026-09-17, see
        //  radiativeRelaxNo() in heterogeneousPyrolysisModel.H)
        scalar maxRad_;
        //- Maximum gas-phase relaxation number (2026-09-21, see
        //  gasRelaxNo() in heterogeneousPyrolysisModel.H)
        scalar maxGas_;

    // Reference to main region fields

        //- InitialPorosity
        volScalarField porosity0_;

        //- Porosity ceiling (2026-09-11): the ash stays solid, so the
        //  porosity must never reach 1 while some ash remains. Computed once
        //  at construction from the initial Yash, rho_ and porosity_ (ash mass
        //  per unit cell volume, constant since the ash is inert). Replaces the
        //  hard-coded clamp 0.9999 -> 1.0 of evolvePorosity().
        volScalarField porosityMax_;

        volScalarField voidFraction_;

        //- Dadiation
        volScalarField& radiation_;

        //- Coefficient of Ts^4 in radiation_ (heterogeneousRadiationModel::
        //  solidRp()) -- used to reintroduce the radiative sink into TEqn
        //  implicitly instead of as a raw explicit source. Reintroduction
        //  fix, 2026-09-14, see the development notes.
        volScalarField& radiationRp_;

    // Fields

        //- Sensible enthalpy gas flux [J/m2/s]
        volScalarField phiHsGas_;

        //- Heat release [J/s/m3]
        volScalarField chemistrySh_;

        //- Porosity growth [1/s]
        volScalarField porositySource_;

        //- Masks for equations
        volScalarField whereIs_;

        volScalarField whereIsNot_;

        volScalarField heatUpGas_;

        //- Reference to heat transport model
        autoPtr<heatTransferModel> HTmodel_;

        volScalarField CONV_;

        volScalarField radiationSh_;

        volTensorField anisotropyK_;

        //- border cells inside solid body
        volScalarField surfF_;

        //- Held igniter (2026-09-10): the cells of igniterCellSet
        //  (pyrolysisProperties) are held at igniterTemperature_ while
        //  time < igniterDuration_, then released. The standard fvOption
        //  fixedTemperatureConstraint cannot be used: it constrains an
        //  enthalpy field through a basicThermo object, not the solid
        //  temperature Ts, and silently has no effect.
        labelList igniterCells_;
        scalar igniterTemperature_;
        scalar igniterDuration_;

    // Checks

        //- Cumulative lost mass of the condensed phase [kg]
        dimensionedScalar lostSolidMass_;

        //- Cumulative mass generation of the gas phase [kg]
        dimensionedScalar addedGasMass_;

        //- Total mass gas flux at the pyrolysing walls [kg/s]
        dimensionedScalar totalGasMassFlux_;

        //- Total heat release rate [J/s]
        dimensionedScalar totalHeatRR_;

        //- scalar chemistryTime
        scalar timeChem_;

    // Bed collapse / shrinkage (2026-09-26, Yang et al. 2005 Fuel 84:2026
    // eq. 7). OPT-IN : pyrolysisCoeffs { bedCollapse on; ... }. With
    // bedCollapse off (default) none of this is used and the model
    // behaves exactly as before (phiS0 not even created).

        //- Master switch (default off)
        Switch bedCollapse_;

        //- Shrinkage coefficient a (h = dz*[phiS0 - a(phiS0-phis)]/phiS0Ref)
        scalar shrinkA_;

        //- Simulated time between two collapses [s]
        scalar collapseInterval_;

        //- Name of the cellZone holding the bed
        word bedZoneName_;

        //- Reference (initial) bed solid fraction
        scalar phiS0Ref_;

        //- Top cell of a column whose solid fraction after remap is below
        //  this value is merged (conservatively) into the cell below
        scalar collapseMinPhis_;

        //- floor(t/collapseInterval) at the last collapse
        scalar lastCollapseIndex_;

        //- Tracer : initial solid volume fraction of the matter now in the
        //  cell (remapped with the matter, never reacts). AUTO_WRITE.
        autoPtr<volScalarField> phiS0Ptr_;

        //- Local cells of the bed zone and their global (column,layer) slot
        labelList bedCells_;
        labelList bedSlot_;

        //- Global slot data (identical on all processors)
        label nBedSlots_;
        labelListList colSlots_;
        scalarField slotZlo_;
        scalarField slotZhi_;
        scalarField slotV_;


    // Protected member functions

        //- Read control parameters
        bool read();

        //- Read control parameters from dict
        bool read(const dictionary& dict);

        //- Mass check
        void calculateMassTransfer();

        // Equations

            //- Solve energy
            void solveEnergy();

            //- Solve solid species mass conservation
            void solveSpeciesMass();

        //- Update whereIs_/whereIsNot_/surfF_ from porosity_ (factored
        //  out of evolvePorosity(), 2026-09-26)
        void updateSolidMasks();

        //- Bed collapse : read controls, build column/layer indexing
        void initBedCollapse();

        //- Bed collapse : conservative 1D remap of the solid along columns
        void collapseBed();

public:

    //- Runtime type information
    TypeName("volPyrolysis");


    // Constructors

        //- Construct from type name, mesh and main region gas list
        volPyrolysis
        (
            const word& modelType,
            const fvMesh& mesh,
            HGSSolidThermo& solidThermo,
            psiReactionThermo& gasThermo,
            volScalarField& whereIs
        );

        //- Construct from type name, mesh, main region gas list and radiation
        volPyrolysis
        (
            const word& modelType,
            const fvMesh& mesh,
            HGSSolidThermo& solidThermo,
            psiReactionThermo& gasThermo,
            volScalarField& whereIs,
            volScalarField& radiation,
            volScalarField& radiationRp
        );

    //- Destructor
    virtual ~volPyrolysis();

    // Member Functions

        // Access

            //- Fields

                //- Return density [kg/m3]
                virtual volScalarField& rhoConst() const;

                //- Return density [kg/m3]
                virtual volScalarField& rho();

                //- Return const temperature [K]
                virtual const volScalarField& T() const;

                //- Return specific heat capacity [J/kg/K]
                virtual const tmp<volScalarField> Cp() const;

                //- Return the region absorptivity [1/m]
                virtual const volScalarField& kappa() const;

                //- Return the region thermal conductivity [W/m/k]
                virtual const volScalarField& K() const;

                //- Return list of cells in surface layer inside solid
                virtual const volScalarField& surf() const;

        // Parametrs from main region

            // Solution parameters

                //- Return the number of non-orthogonal correctors
                inline label nNonOrthCorr() const;

                //- Return max diffusivity allowed in the solid
                virtual scalar maxDiff() const;

        // Helper functions

            //- Mean diffusion number of the solid region
            virtual scalar solidRegionDiffNo() const;

            //- Radiative relaxation number (2026-09-17, see
            //  heterogeneousPyrolysisModel.H)
            virtual scalar radiativeRelaxNo() const;

            //- Max radiative relaxation number allowed
            virtual scalar maxRadDiff() const;

            //- Gas-phase relaxation number (2026-09-21, see
            //  heterogeneousPyrolysisModel.H)
            virtual scalar gasRelaxNo();

            //- Max gas-phase relaxation number allowed
            virtual scalar maxGasDiff() const;

            //- Returm time modification factor
            virtual scalar maxTime() const;

            //- Returm time modification factor
            virtual Switch equilibrium() const;

        // Source fields (read/write access)

       // Evolution

            //- Pre-evolve region
            virtual void preEvolveRegion();

            //- Evolve the pyrolysis equations
            virtual void evolveRegion();

            //- Re-solve the energy equation only (radiation<->Ts
            //  sub-iteration, see heterogeneousPyrolysisModel.H)
            virtual void solveEnergyOnly();

            //- Evolve porosity field
            virtual void evolvePorosity();

        //Comunication with main region

            //- Return total mass source - Eulerian phase only
            virtual Foam::tmp<Foam::volScalarField> Srho() const;

            //- Return mass source for specie i - Eulerian phase only
            virtual Foam::tmp<Foam::volScalarField> Srho
            (
                const label i
            ) const;

            //- Return enthalpy source - Eulerian phase only
            virtual Foam::tmp<Foam::volScalarField> heatTransfer();

            //- Return CONV
            virtual Foam::tmp<Foam::volScalarField> CONV() const;

            //- Return enthalpy needed to heat up gas
            virtual Foam::tmp<Foam::volScalarField> heatUpGasCalc() const;

            //- Return enthalpy needed to heat up gas
            virtual Foam::tmp<Foam::volScalarField> heatUpGas() const;

            //- Return enthalpy needed to heat up gas  GK  or to pyrolyse the solid into gas ???
            virtual Foam::tmp<Foam::volScalarField> solidChemistrySh() const;

        // I-O

            //- Provide some feedback
            virtual void info() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heterogeneousPyrolysisModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volPyrolysisI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
