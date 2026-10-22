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

// NaN check (2026-10-13). A NaN can spread to ~1100 cells in one step, and a check after
// EEqn alone does not show where it starts. Each stage of the time step is checked (fields,
// and matrices before the solution); the run stops at the FIRST stage where a non-finite
// value appears and prints the first cells concerned with their context. Finite values
// outside the physical range are only reported (warning, bounded number of lines).
// Each check performs a reduce: all processors must call the same checks in the same order.

label dbgWarnLines = 0;
const label dbgWarnMax = 400;

auto dbgCtx = [&](const label c)
{
    const volScalarField& TsDbg = pyrolysisZone.T();
    Pout<< "DBG2     ctx C " << mesh.C()[c] << " p " << p[c] << " rho " << rho[c]
        << " T " << T[c] << " Told " << T.oldTime()[c] << " Ts " << TsDbg[c]
        << " eps " << porosityF[c] << " U " << U[c] << " Srho " << Srho[c]
        << " heatTransfer " << heatTransfer[c] << nl << "DBG2     Y:";
    forAll(Y, i) { Pout<< " " << Y[i].name() << "=" << Y[i][c]; }
    Pout<< nl << "DBG2     Yold:";
    forAll(Y, i) { Pout<< " " << Y[i].name() << "=" << Y[i].oldTime()[c]; }
    Pout<< endl;
};

// Returns the global number of non-finite cells; also reports those outside [lo, hi].
auto dbgField = [&](const word& stage, const volScalarField& f, const scalar lo, const scalar hi) -> label
{
    label nNan = 0;
    label nOut = 0;
    forAll(f, c)
    {
        const scalar v = f[c];
        if (!std::isfinite(v))
        {
            if (nNan < 5)
            {
                Pout<< "DBG2 NAN stage " << stage << " field " << f.name() << " cell " << c
                    << " val " << v << endl;
                dbgCtx(c);
            }
            nNan++;
        }
        else if ((v < lo || v > hi) && dbgWarnLines < dbgWarnMax && nOut < 3)
        {
            Pout<< "DBG2 WARN t " << runTime.timeName() << " stage " << stage << " field " << f.name()
                << " cell " << c << " val " << v << endl;
            dbgCtx(c);
            nOut++;
            dbgWarnLines++;
        }
    }
    reduce(nNan, sumOp<label>());
    return nNan;
};

// Matrix before the solution: diagonal, source and off-diagonal coefficients.
auto dbgMatrix = [&](const word& stage, const fvScalarMatrix& M) -> label
{
    label nBad = 0;
    const scalarField& d = M.diag();
    const scalarField& s = M.source();
    forAll(d, c)
    {
        if (!std::isfinite(d[c]) || !std::isfinite(s[c]) || d[c] <= 0)
        {
            if (nBad < 5)
            {
                Pout<< "DBG2 MATRIX stage " << stage << " cell " << c << " diag " << d[c]
                    << " source " << s[c] << " V " << mesh.V()[c] << endl;
                dbgCtx(c);
            }
            nBad++;
        }
    }
    if (M.hasUpper())
    {
        const scalarField& u = M.upper();
        forAll(u, facei)
        {
            if (!std::isfinite(u[facei]))
            {
                if (nBad < 5)
                {
                    Pout<< "DBG2 MATRIX stage " << stage << " upper face " << facei
                        << " owner " << mesh.owner()[facei] << " neigh " << mesh.neighbour()[facei]
                        << " val " << u[facei] << endl;
                    dbgCtx(mesh.owner()[facei]);
                }
                nBad++;
            }
        }
    }
    if (M.hasLower())
    {
        const scalarField& l = M.lower();
        forAll(l, facei)
        {
            if (!std::isfinite(l[facei]))
            {
                if (nBad < 5)
                {
                    Pout<< "DBG2 MATRIX stage " << stage << " lower face " << facei
                        << " owner " << mesh.owner()[facei] << " val " << l[facei] << endl;
                    dbgCtx(mesh.owner()[facei]);
                }
                nBad++;
            }
        }
    }
    reduce(nBad, sumOp<label>());
    return nBad;
};

auto dbgStop = [&](const word& stage, const label n)
{
    if (n > 0)
    {
        Info<< "DBG2 FIRST BAD stage " << stage << " nBad " << n << " t " << runTime.timeName() << endl;
        runTime.writeNow();
        FatalErrorInFunction << "DBG2 stop at stage " << stage << exit(FatalError);
    }
};
