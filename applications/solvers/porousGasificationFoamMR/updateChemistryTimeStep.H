if (adjustTimeStep and solidChemistryTimeStepControl)
{

    const scalar TsolidChemistry = pyrolysisZone.maxTime();
    const scalar dt0 = runTime.deltaTValue();

    if (min(TsolidChemistry,dt0) == TsolidChemistry)
    {
        Info << "timestep corrected for heterogeneous reactions to: " << TsolidChemistry << "[s]" << endl;
    }

    runTime.setDeltaT
    (
        min
        (
            TsolidChemistry,
            dt0
        )
    );

}

// Gas-phase reactive timestep control. Previously gasChemistryTimeStepControl
// was read (readChemistryTimeControls.H) but never consulted anywhere in the
// non-LTS (transient) branch -- deltaT was sized purely from flow Courant +
// solid diffusion + solid chemistry, leaving the gas Arrhenius ODE (seulex)
// to receive an outer dt with no regard for its own stiffness. Root-caused
// via SIGFPE in IrreversibleReaction::kf() on the very first outer step
// (2026-09-10, see the development notes): kf() = A*exp(-Ta/T) diverges
// when seulex's internal adaptive substepping, integrating over a dt far
// larger than the gas reaction's natural timescale, produces a non-physical
// trial state. Confirmed independent of h/SAV and of the region-coupling
// gas preheat (identical crash with both removed).
if (adjustTimeStep and gasChemistryTimeStepControl)
{
    const dictionary& chemDict =
        mesh.lookupObject<dictionary>("chemistryProperties");

    scalar TgasChemistry = great;

    if (runTime.timeIndex() == runTime.startTimeIndex())
    {
        // reaction->Qdot() is exactly zero before gas chemistry has ever
        // been solved once, so a heat-release-based bound (below) cannot
        // protect this very first reactive step -- it would see Qdot=0
        // everywhere and return no restriction at all, which is exactly
        // the case that crashed. Bound this one step directly from the
        // same seed value basicChemistryModel/basicSolidChemistryModel
        // already read for their own internal deltaTChem_ initialisation
        // (deltaTChemIni_ <- "initialChemicalTimeStep").
        TgasChemistry = chemDict.lookupOrDefault<scalar>
        (
            "initialChemicalTimeStep",
            great
        );
    }
    else
    {
        // Heat-release-rate timescale: bound dt so no cell's temperature
        // would change by more than a fraction gasAlphaTemp of itself in
        // one step, from the *previous* step's Qdot -- same criterion as
        // reactingFoam's LTS setRDeltaT.H, reduced here to a single global
        // scalar (gMax over all cells) to match this solver's non-LTS,
        // single global deltaT architecture.
        scalar alphaTemp
        (
            chemDict.lookupOrDefault<scalar>("gasAlphaTemp", 0.05)
        );

        volScalarField::Internal rDeltaTT
        (
            mag(reaction->Qdot())/(alphaTemp*rho*thermo.Cp()*T)
        );

        TgasChemistry = 1.0/(gMax(rDeltaTT.field()) + vSmall);
    }

    const scalar dt0 = runTime.deltaTValue();

    if (min(TgasChemistry,dt0) == TgasChemistry)
    {
        Info << "timestep corrected for gas-phase chemistry to: " << TgasChemistry << "[s]" << endl;
    }

    runTime.setDeltaT
    (
        min
        (
            TgasChemistry,
            dt0
        )
    );

}
