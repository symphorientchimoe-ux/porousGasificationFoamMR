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
    Local particle diameter of a bed of shrinking spheres (article 2, 2026-10-30).
    The number of particles in a cell does not change during conversion (no
    fragmentation). It is carried by the solid-volume tracer phiS0 of the bed
    collapse (initial solid volume fraction, remapped with the particles), or
    by 1 - initialPorosity when there is no collapse. Hence
        d = d0 * cbrt( (1 - porosity) / phiS0 ),   clamped to [dMin, d0],
    and the specific surface of the bed is a_v = 6 (1 - porosity) / d.
    In a cell with no solid (phiS0 = 0) the surface is zero.

\*---------------------------------------------------------------------------*/

#ifndef shrinkingParticle_H
#define shrinkingParticle_H

#include "volFields.H"

namespace Foam
{
namespace shrinkingParticle
{

//- Initial solid volume fraction of cell celli (phiS0 if registered, else 1 - initialPorosity)
inline scalar phiS0
(
    const fvMesh& mesh,
    const volScalarField& initialPorosity,
    const label celli
)
{
    if (mesh.foundObject<volScalarField>("phiS0"))
    {
        return mesh.lookupObject<volScalarField>("phiS0")[celli];
    }
    return 1.0 - initialPorosity[celli];
}

//- Particle diameter [m]
inline scalar diameter
(
    const scalar d0,
    const scalar dMin,
    const scalar solidFraction,
    const scalar phiS0
)
{
    if (phiS0 <= small || solidFraction <= 0)
    {
        return dMin;
    }
    return min(d0, max(dMin, d0*cbrt(min(solidFraction/phiS0, 1.0))));
}

//- Specific surface of the bed [1/m]
inline scalar specificSurface(const scalar solidFraction, const scalar d)
{
    return solidFraction > 0 ? 6.0*solidFraction/d : 0.0;
}

} // End namespace shrinkingParticle
} // End namespace Foam

#endif
