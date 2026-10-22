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
    Foam::heterogeneousRadiationModel

Description
    Top level model for heterogeneous radiation modelling

SourceFiles
    heterogeneousRadiationModel.C

\*---------------------------------------------------------------------------*/

#ifndef heterogeneousRadiationModel_H
#define heterogeneousRadiationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "fvMatricesFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class basicThermo;
class fvMesh;

namespace radiationModels
{
class heterogeneousAbsorptionEmissionModel;
class scatterModel;
class sootModel;
}

/*---------------------------------------------------------------------------*\
                 Class heterogeneousRadiationModel Declaration
\*---------------------------------------------------------------------------*/

class heterogeneousRadiationModel
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;

        //- Reference to the gas temperature field
        const volScalarField& T_;

        //- Reference to the solid temperature field
        const volScalarField& Ts_;

        //- Radiation model dictionary
        dictionary coeffs_;

        //- Radiation solver frequency - number flow solver iterations per
        //  radiation solver iteration
        label solverFreq_;

        //- Flag to enable radiation model to be evaluated on first iteration
        bool firstIter_;


        // References to the radiation sub-models

            //- Absorption/emission model
            autoPtr<radiationModels::heterogeneousAbsorptionEmissionModel>
                heterogeneousAbsorptionEmission_;

            //- Scatter model
            autoPtr<radiationModels::scatterModel> scatter_;

            //- Soot model
            autoPtr<radiationModels::sootModel> soot_;


private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Initialise
        void initialise();


public:

    //- Runtime type information
    TypeName("heterogeneousRadiationModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            heterogeneousRadiationModel,
            T,
            (
                const volScalarField& T
            ),
            (T)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            heterogeneousRadiationModel,
            porosity,
            ( 
                const volScalarField& T,
                const volScalarField& porosityF,
                const volScalarField& surfF,
                const volScalarField& Ts
            ),
            (T,porosityF,surfF,Ts)
        );

    // Constructors

        //- Null constructor
        heterogeneousRadiationModel(const volScalarField& T);

        //- Construct from components
        heterogeneousRadiationModel(const word& type, const volScalarField& T);

        //- Construct from components with a distinct solid temperature
        //  (2026-09-17, see the development notes) : the 2-arg
        //  constructor above sets Ts_(T) -- fine for the null/no-radiation
        //  path, but heterogeneousP1 was (wrongly) using it too, silently
        //  binding Ts_ to the GAS temperature instead of the solid Ts
        //  parameter it actually receives. Root cause of the whole
        //  borderAs investigation (runs 8-26) : solidRadiation was always
        //  computed from Tgas, never from the real (possibly crashed,
        //  possibly very negative) solid temperature -- no negative
        //  feedback from the solid's own state, so nothing bounded its
        //  radiative sink once it started cooling.
        heterogeneousRadiationModel
        (
            const word& type,
            const volScalarField& T,
            const volScalarField& Ts
        );

        //- Construct from porosity, surface and solid temperature
        heterogeneousRadiationModel
        (
            const volScalarField& T,
            const volScalarField& porosityF,
            const List<label>& surfF,
            const volScalarField& Ts
        );

    // Selectors

        //- Return a reference to the selected radiation model
        static autoPtr<heterogeneousRadiationModel> New(const volScalarField& T);

        //- Return a reference to the selected radiation model
        static autoPtr<heterogeneousRadiationModel> New
        (
            const volScalarField& T,
            const volScalarField& porosityF,
            const volScalarField& surfF,
            const volScalarField& Ts
        );


    //- Destructor
    virtual ~heterogeneousRadiationModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();

            //- Solve radiation equation(s)
            virtual void calculate() = 0;

            //- Read radiationProperties dictionary
            virtual bool read() = 0;


        // Access

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const = 0;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const = 0;

            //- Energy source term
            virtual tmp<fvScalarMatrix> Sh
            (
                const basicThermo& thermo,
                const volScalarField& he
            ) const;

            //- Temperature source term
            virtual tmp<fvScalarMatrix> ST
            (
                const dimensionedScalar& rhoCp,
                volScalarField& T
            ) const;

            //- Enthalpy source term for solid
            virtual tmp<volScalarField> solidSh() const = 0;

            //- Coefficient of Ts^4 in solidSh() (>=0), for implicit
            //  linearisation of the radiative sink in the solid energy
            //  equation: solidSh(Ts) == solidRu - solidRp()*Ts^4, with
            //  solidRu independent of Ts. Reintroduced 2026-09-14 as the
            //  real fix for the heterogeneousP1 divergence (see
            //  the development notes) -- solidSh() alone was always
            //  applied as a fully explicit source in
            //  volPyrolysis::solveEnergy(), so nothing prevented a large
            //  sigma*Ts^4 sink from overshooting T below 0K in a single
            //  timestep once combustion reached high temperature.
            virtual tmp<volScalarField> solidRp() const = 0;

            //- Access to absorptionEmission model
            const radiationModels::heterogeneousAbsorptionEmissionModel&
                heterogeneousAbsorptionEmission() const;

            //- Access to soot Model
            const radiationModels::sootModel& soot() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const heterogeneousRadiationModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
