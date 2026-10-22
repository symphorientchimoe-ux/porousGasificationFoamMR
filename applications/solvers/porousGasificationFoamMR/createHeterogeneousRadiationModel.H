    autoPtr<heterogeneousRadiationModel> radiation
    (
        heterogeneousRadiationModel::New(thermo.T(),porosityF,pyrolysisZone.surf(),pyrolysisZone.T())
    );

    // Number of radiation<->Ts sub-iterations per timestep (2026-09-16, see
    // the development notes) : default 1 keeps the pre-existing
    // behaviour (radiation.H called once, before pyrolysisZone.evolve(),
    // Ts sees a G frozen from the previous timestep) unchanged unless
    // explicitly raised in radiationProperties. >1 re-runs radiation-
    // >correct() + pyrolysisZone.solveEnergyOnly() with the just-updated
    // Ts_, tightening the coupling during fast transients (e.g. ignition).
    // Read from the already-constructed `radiation` object itself
    // (heterogeneousRadiationModel publicly inherits IOdictionary, already
    // MUST_READ and registered under the name "radiationProperties") --
    // NOT a fresh IOdictionary of the same name, which would collide with
    // that existing registration (regIOobject: object already registered).
    label nRadCorrectors(radiation().lookupOrDefault<label>("nRadCorrectors", 1));
