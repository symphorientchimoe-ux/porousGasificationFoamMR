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
    Gas-solid mass transfer of a bed of shrinking spheres (article 2, 2026-10-30).
    ST = h_m a_v with a_v = 6 (1 - porosity)/d and d from shrinkingParticle.H.
    hModel constant: h_m read from the dictionary (only the surface changes);
    hModel WakaoFunazkri: Sh = 2 + 1.1 Re^0.6 Sc^(1/3), Re = rho |U| d/mu,
    Sc = mu/(rho D), D = Dref (T/Tref)^1.75, h_m = Sh D/d.

\*---------------------------------------------------------------------------*/

#ifndef shrinkingSphereST_H
#define shrinkingSphereST_H

#include "specieTransferModel.H"

namespace Foam
{

class shrinkingSphereST
:
    public specieTransferModel
{
        scalar d0_;
        scalar dMin_;
        Switch correlated_;
        scalar h_;
        scalar Dref_;
        scalar Tref_;

public:

    TypeName("shrinkingSphereST");

    shrinkingSphereST(const volScalarField& por, const volScalarField& por0);

    //- Selector used by the run-time selection table (declareRunTimeNewSelectionTable)
    static autoPtr<shrinkingSphereST> New
    (
        const volScalarField& por,
        const volScalarField& por0
    );

    virtual ~shrinkingSphereST()
    {}

    //- Return the gas-solid mass transfer field h_m a_v [1/s]
    virtual tmp<volScalarField> ST() const;

    virtual bool read();
};

} // End namespace Foam

#endif
