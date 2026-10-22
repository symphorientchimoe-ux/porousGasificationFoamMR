// 2026-10-12: RRg, hence Srho, is already per unit total volume (Ym = Ys*rhos*(1-eps)); the factor (1-porosityF)
// removed eps*Srho (41 % of the released gas, paper Section 3.7). Consistent with YEqn (full Srho(i)).
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Global
    rhoEqn

Description
    Solve the continuity for density.

\*---------------------------------------------------------------------------*/
// eqZx2uHGn042

Srho = pyrolysisZone.Srho();

{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(porosityF, rho)
      + fvc::div(phi)
      ==
        Srho   // 2026-10-12: no (1 - porosityF) factor, see the header
    );

    fvOptions.constrain(rhoEqn);

    rhoEqn.solve();

    fvOptions.correct(rho);
}

// ************************************************************************* //
