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
Class
    Foam::combustionModels::strainMixingLimited

Description
    Mixing-limited combustion model for laminar (no RAS/LES) flows.

    Same structure as the standard PaSR combustion model
    (kappa = tc/(tc+tmix), R = kappa*laminar::R, Qdot = kappa*laminar::Qdot),
    but PaSR's turbulence-based mixing time
        tmix_PaSR = Cmix*sqrt(muEff/rho/epsilon)
    requires a RAS/LES epsilon field that this case does not have
    (simulationType laminar; epsilon=0 in a laminar momentumTransportModel,
    which would make tmix -> infinity and kappa -> 0, i.e. combustion fully
    extinguished -- not usable as-is).

    Here the mixing time is instead estimated directly from the RESOLVED
    velocity field's strain rate, available in any laminar simulation:

        strainRate = mag(symm(fvc::grad(U)))                      [1/s]
        tmix        = Cmix / (strainRate + smallRate)              [s]
        kappa       = tc / (tc + tmix)

    This is grounded in flamelet-type reasoning (the local strain rate
    controls the scalar dissipation rate, hence the rate at which a
    fuel/oxidiser interface can mix, in both laminar and turbulent
    diffusion flames) rather than in an empirically-tuned rate constant.
    Cmix is a single O(1) model constant, taken from the literature range
    used for the standard PaSR model (not fitted to match a target
    temperature) -- see constant/fluidZone/combustionProperties.

    U is looked up from the mesh object registry (mesh().lookupObject),
    since the combustionModel base class only exposes rho()/phi() from the
    turbulence model, not U directly.

SourceFiles
    strainMixingLimited.C

\*---------------------------------------------------------------------------*/

#ifndef strainMixingLimited_H
#define strainMixingLimited_H

#include "laminar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                   Class strainMixingLimited Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class strainMixingLimited
:
    public laminar<ReactionThermo>
{
    // Private Data

        //- Mixing constant (O(1), literature PaSR range -- not tuned)
        scalar Cmix_;

        //- Mixing parameter kappa = tc/(tc+tmix), in [0,1]
        volScalarField kappa_;

        //- Mixing time scale 1/strainRate, written for diagnostics
        volScalarField tmix_;


public:

    //- Runtime type information
    TypeName("strainMixingLimited");


    // Constructors

        //- Construct from components
        strainMixingLimited
        (
            const word& modelType,
            const ReactionThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        strainMixingLimited(const strainMixingLimited&);


    //- Destructor
    virtual ~strainMixingLimited();


    // Member Functions

        //- Correct combustion rate
        virtual void correct();

        //- Fuel consumption rate matrix.
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s^3]
        virtual tmp<volScalarField> Qdot() const;

        //- Update properties from given dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const strainMixingLimited&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "strainMixingLimited.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
