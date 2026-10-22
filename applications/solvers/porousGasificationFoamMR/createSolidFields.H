// Refractory thermal fields, region solidZone (2026-09-09). Pure conduction
// (no species, no flow), as in chtMultiRegionFoam. Uses the standard OpenFOAM
// solidThermo (constant/solidZone/thermophysicalProperties, heSolidThermo/
// pureMixture/constIso/eConst/rhoConst): the interface condition
// compressible::turbulentTemperatureCoupledBaffleMixed (EEqn.H on the fluid
// side, T on the solid side) finds kappa() through temperatureCoupledBase,
// which requires a registered solidThermo object, not a free kappa field.

Info<< "*** Reading solid mesh thermophysical properties for region "
    << solidMesh.name() << nl << endl;

// Foam::solidThermo qualified explicitly: createPyrolysisModel.H (included
// earlier, unchanged from porousGasificationFoam) already declares a local
// reference named "solidThermo" (HGSSolidThermo&), which shadows the class
// name for unqualified lookup.
autoPtr<Foam::solidThermo> solidThermoPtr(Foam::solidThermo::New(solidMesh));
Foam::solidThermo& thermoSolid = solidThermoPtr();

volScalarField& eSolid = thermoSolid.he();
volScalarField& rhoSolid = thermoSolid.rho();

// ************************************************************************* //
