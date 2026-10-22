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
    Foam::combustionModels::laminarEDM

Description
    Finite-rate / Eddy-Dissipation combustion model (Magnussen & Hjertager
    1977) for laminar (no RAS/LES) flows. Added 2026-09-25 (test56).

    For each gas reaction r (fuel f_r + s_r kg O2 per kg fuel -> products):

        R_kin,r = instantaneous Arrhenius fuel consumption rate
                  (chemistry calculateRR, same reactions/A/Ta as laminar)
        R_EDM,r = A*rho*strainRate*min(Y_f, Y_O2/s_r, B*Y_P,r/(1 + s_r))
        R_r     = min(R_kin,r, R_EDM,r)                    [kg fuel/m3/s]

    with A=4, B=0.5 (Magnussen & Hjertager 1977, standard values, not
    tuned). The classical EDM mixing rate epsilon/k needs a RAS model this
    case does not have (simulationType laminar; kEpsilon was ~1000x slower,
    test25), so it is replaced by the resolved strain rate
        strainRate = sqrt(2)*mag(symm(grad(U)))            [1/s]
    -- same substitution as strainMixingLimited, but with the true EDM
    min() form instead of the PaSR kappa=tc/(tc+tmix) blend.

    The gas-phase ODE integration of laminar::correct() is NOT called: the
    rates are the instantaneous finite-rate/EDM rates (as in the classical
    finite-rate/eddy-dissipation formulation), so the Rosenbrock gas ODE
    solver is no longer used.

    Numerical safeguard (not physics): each R_r is also capped so that a
    single time step cannot consume more fuel or O2 than present in the
    cell (realizability, rho*min(Y_f, Y_O2/s_r)/deltaT, and a common scaling
    when the reactions together would exceed the available O2). The
    fraction of cells where this cap is active is printed every step.

    Per-reaction data, in the SAME ORDER as the reactions in
    constant/chemistryProperties, from combustionProperties:

        laminarEDMCoeffs
        {
            A        4;
            B        0.5;
            oxidant  O2;
            reactions
            (
                COoxidation  { fuel CO;  s 0.5714; products (CO2); }
                ...
            );
        }

SourceFiles
    laminarEDM.C

\*---------------------------------------------------------------------------*/

#ifndef laminarEDM_H
#define laminarEDM_H

#include "laminar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                         Class laminarEDM Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class laminarEDM
:
    public laminar<ReactionThermo>
{
    // Private Data

        //- Magnussen-Hjertager constants
        scalar A_;
        scalar B_;

        //- Oxidant specie index
        label oxidanti_;

        //- Per reaction: fuel specie index
        labelList fueli_;

        //- Per reaction: stoichiometric O2/fuel mass ratio
        scalarList s_;

        //- Per reaction: product specie indices
        List<labelList> productsi_;

        //- Specie mass reaction rates [kg/m3/s]
        PtrList<volScalarField::Internal> RR_;

        //- Per reaction: R_r/R_kin,r (1 = kinetics-limited), diagnostics
        PtrList<volScalarField> kappa_;


    // Private Member Functions

        //- Read the per-reaction data
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("laminarEDM");


    // Constructors

        //- Construct from components
        laminarEDM
        (
            const word& modelType,
            const ReactionThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        laminarEDM(const laminarEDM&);


    //- Destructor
    virtual ~laminarEDM();


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
        void operator=(const laminarEDM&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "laminarEDM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
