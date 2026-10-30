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
Description
    Gas-solid heat exchange of a bed of shrinking spheres (article 2, 2026-10-30).
    CONV = h a_v with a_v = 6 (1 - porosity)/d and d from shrinkingParticle.H.
    hModel constant: h read from the dictionary (only the surface changes);
    hModel WakaoKaguei: Nu = 2 + 1.1 Re^0.6 Pr^(1/3), Re = rho |U| d/mu,
    h = Nu kappa/d with kappa = Cp alpha (thermo:alpha = kappa/Cp in OpenFOAM 8).
    With d0 and the initial porosity of constCONV, the initial a_v equals its SAV.

\*---------------------------------------------------------------------------*/

#ifndef shrinkingSphereCONV_H
#define shrinkingSphereCONV_H

#include "heatTransferModel.H"
#include "fluidThermo.H"

namespace Foam
{

class shrinkingSphereCONV
:
    public heatTransferModel
{
    // Private data

        //- Initial and minimum particle diameters [m]
        scalar d0_;
        scalar dMin_;

        //- Correlated h (WakaoKaguei) or constant h
        Switch correlated_;

        //- Constant heat transfer coefficient [W/m2/K]
        scalar h_;

        const volVectorField& Up_;
        const volScalarField& rhop_;
        const volScalarField& alphap_;
        const volScalarField& mup_;
        const fluidThermo& thermop_;

public:

    TypeName("shrinkingSphereCONV");

    shrinkingSphereCONV
    (
        const volScalarField& porosity,
        const volScalarField& initialPorosity
    );

    //- Selector used by the run-time selection table (declareRunTimeNewSelectionTable)
    static autoPtr<shrinkingSphereCONV> New
    (
        const volScalarField& porosity,
        const volScalarField& initialPorosity
    );

    virtual ~shrinkingSphereCONV()
    {}

    //- Return the field for heat transfer between solid/gas [W/m3/K]
    virtual tmp<volScalarField> CONV() const;

    //- Read heatTransferProperties dictionary
    virtual bool read();
};

} // End namespace Foam

#endif
