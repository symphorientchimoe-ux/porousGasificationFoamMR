// 2026-10-12: RRg, hence Srho, is already per unit total volume (Ym = Ys*rhos*(1-eps)); the factor (1-porosityF)
// removed eps*Srho (41 % of the released gas, paper Section 3.7). Consistent with YEqn (full Srho(i)).
rho = thermo.rho();
dbgStop("G_pEqn_thermoRho", dbgField("G_pEqn_thermoRho", rho, 0.02, 3.0));   // NaN check

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pimple.nCorrPiso() <= 1)
{
    tUEqn.clear();
}

if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            fvc::flux(HbyA)
          + MRF.zeroFilter
            (
                rhorAUf*fvc::ddtCorr(rho, U, phi)/fvc::interpolate(rho)
            )
        )
    );

    MRF.makeRelative(fvc::interpolate(psi), phid);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rhorAUf, p)
         ==
            fvOptions(psi, p, rho.name())
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rho*HbyA)
        )
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(porosityF * psi, p)
          + fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p)
         ==
            fvOptions(psi * porosityF, p, rho.name())
          + Srho   // 2026-10-12: no (1 - porosityF) factor, see the header
        );

        // v1.1 (2026-10-26, paper Section 3.3): ddt(porosityF*psi, p) acts on a
        // temporary whose old-time value is its current value, so it gives
        // eps psi dp/dt only and loses p d(eps psi)/dt (porosity change and
        // thermal expansion). The missing part, explicit in p.oldTime(), makes
        // the storage term (eps psi p - (eps psi p)^o)/dt, as in rhoEqn.
        if (conservativeDdt)
        {
            pEqn += fvc::ddt(porosityF, psi)*p.oldTime();
        }

        dbgStop("G_pEqn_matrix", dbgMatrix("G_pEqn_matrix", pEqn));   // NaN check

        pEqn.solve();

        dbgStop("G_pEqn_p", dbgField("G_pEqn_p", p, 9e4, 1.1e5));   // NaN check

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
dbgStop("G_pEqn_rhoEqn", dbgField("G_pEqn_rhoEqn", rho, 0.02, 3.0));   // NaN check
#include "compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
if (!pimple.finalPimpleIter())
{
    p.relax();
    rho.relax();
}

Info<< "rho max/min : " << max(rho).value()
    << " " << min(rho).value() << endl;

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

// H5 (2026-09-26, code analysis): K is created once in createFields.H and was only
// updated in pcEqn.H (unused), so fvc::ddt(...,K) was 0 and div(phi,K) used the K of
// the restart. Opt-in switch so running/restarted cases are unchanged:
// system/controlDict "updateK yes;" (v1.1, 2026-10-26: default yes, updateK no = v1.0)
if (runTime.controlDict().lookupOrDefault<Switch>("updateK", true))
{
    K = 0.5*magSqr(U);
}

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
