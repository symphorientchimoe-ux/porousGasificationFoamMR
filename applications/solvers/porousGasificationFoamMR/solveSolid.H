// Conduction in solidZone (refractory) with the standard solidThermo (see
// createSolidFields.H). The fluidZone/solidZone interface is coupled by the
// boundary condition compressible::turbulentTemperatureCoupledBaffleMixed
// (or its radiative variant) on T on both sides (0/fluidZone/T, 0/solidZone/T).
{
    fvScalarMatrix eSolidEqn
    (
        fvm::ddt(rhoSolid, eSolid)
      - fvm::laplacian(thermoSolid.alpha(), eSolid)
    );

    eSolidEqn.relax();
    eSolidEqn.solve();
}

thermoSolid.correct();

Info<< " T (solidZone) min/max = " << gMin(thermoSolid.T()) << ", "
    << gMax(thermoSolid.T()) << endl;

// ************************************************************************* //
