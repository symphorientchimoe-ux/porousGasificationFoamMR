// Two fixed meshes with hard-coded names (2026-09-09); constant/regionProperties
// is not read, the solver has exactly two regions. fluidZone carries the
// porousGasificationFoam physics (plenum, biomass_bed, gas_zone); solidZone is
// the refractory, solved separately by conduction (see createSolidFields.H).

const word fluidZoneName("fluidZone");
const word solidZoneName("solidZone");

Info<< "Create fluid mesh (" << fluidZoneName << ") for time = "
    << runTime.timeName() << nl << endl;

fvMesh mesh
(
    IOobject
    (
        fluidZoneName,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);

Info<< "Create solid mesh (" << solidZoneName << ") for time = "
    << runTime.timeName() << nl << endl;

fvMesh solidMesh
(
    IOobject
    (
        solidZoneName,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);

// ************************************************************************* //
