scalar maxDi = pyrolysis->maxDiff();

// 2026-09-17 (see the development notes) : caps the radiative
// relaxation number (radiativeRelaxNo(), controlDict "maxRa", default 0.3
// if absent) alongside the existing diffusion number -- see
// solidRegionDiffusionNo.H/setMultiRegionDeltaT.H.
scalar maxRa = pyrolysis->maxRadDiff();

// 2026-09-21 (see gasRelaxNo() in heterogeneousPyrolysisModel.H) : caps
// the new gas-phase relaxation number, controlDict "maxGa", default 0.3
// if absent -- same role as maxRa above but protects the explicit
// heatTransfer/heatUpGas source terms in EEqn.H, which never had a dt
// criterion even though the gas thermal mass is far smaller than the
// solid's.
scalar maxGa = pyrolysis->maxGasDiff();

// 2026-09-21 : caps the new velocity-growth number (UaNum, see
// solidRegionDiffusionNo.H) -- same role as maxRa/maxGa above but for the
// momentum equation (UEqn.H/pEqn.H). Instrumentation (tests 11-13) showed
// the JANAF-ceiling divergence on this mesh happens entirely in gas_zone
// (porosityF=1 -- no Darcy resistance there, by design), driven by
// advection (embedded in UEqn's H-operator, feeding HbyA in pEqn.H) that
// keeps amplifying U from one timestep to the next with nothing to check
// it -- not buoyancy (rho*g stayed small, wrong direction to explain the
// growing component) and not the pressure-gradient correction (small
// compared to HbyA). No existing criterion bounds dt by how fast U itself
// is changing (Courant only measures U*dt/dx, not dU/dt). Purely a
// generic flow-field quantity (no pyrolysisModel involved), read directly
// here like maxCo/maxDi. Not yet calibrated.
scalar maxUa = runTime.controlDict().lookupOrDefault<scalar>("maxUa", 0.3);
