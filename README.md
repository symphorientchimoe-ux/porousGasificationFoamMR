# porousGasificationFoamMR

An OpenFOAM-8 solver for the thermochemical conversion of a fixed bed of solid fuel (drying, devolatilisation, char oxidation and gasification, gas-phase combustion) inside a refractory-lined reactor. It is derived from porousGasificationFoam (P.J. Żuk et al., Comput. Phys. Commun. 278 (2022) 108407, https://github.com/btuznik/porousGasificationFoam) and keeps its porous-medium description of the bed.

What this version adds:

- a second mesh region (`solidZone`, the refractory wall) solved by conduction and coupled to the gas region (`fluidZone`) by conjugate temperature and Marshak radiative boundary conditions;
- a laminar eddy-dissipation combustion model, `laminarEDM`, in which the resolved strain rate replaces the turbulent dissipation frequency;
- an Ergun bed resistance recomputed at every time step from the local porosity (viscous term, and optionally the inertial term);
- an opt-in, mass-conserving bed-collapse remap (Yang's shrinkage law);
- extra time-step criteria (refractory diffusion number, radiative, gas-exchange, velocity-growth and gas heat-release limits);
- corrections that close the global mass and energy balances, described in the accompanying paper.

The physics, the corrections and the validation against the counter-current wood-pellet reactor of Kim et al. (Waste Management 49, 2016) are described in the accompanying article (Computer Physics Communications, submitted).

## Requirements

- OpenFOAM-8 (The OpenFOAM Foundation). Version used for the paper: build 8-a86b07b9dabb.
- A C++ compiler supported by OpenFOAM-8 (tested with GCC on Linux x86-64).
- Python 3 with numpy and matplotlib, only for the post-processing scripts of the test cases.

## Building

```sh
source /path/to/OpenFOAM-8/etc/bashrc
./Allwmake
```

The libraries are installed in `$FOAM_USER_LIBBIN` and the solver in `$FOAM_USER_APPBIN`, the standard OpenFOAM user locations. `./Allwclean` removes the build files.

The libraries keep the names of porousGasificationFoam (`libHGSradiationModels`, `libHGSheterogeneousPyrolysisModels`, `libfieldPorosityModel`, `libHGSsolid`, `libHGSSolidThermo`, `libsolidChemistryModel`). If the original porousGasificationFoam is installed in the same `$FOAM_USER_LIBBIN`, the two will overwrite each other: install them under different user directories (set `WM_PROJECT_USER_DIR`, or `FOAM_USER_LIBBIN` and `FOAM_USER_APPBIN`, before building).

## Directory layout

```
applications/solvers/porousGasificationFoamMR   solver (two regions, laminarEDM, Ergun, time-step control)
src/porousGasificationMedia                     libraries: pyrolysis (volPyrolysis, bed collapse), radiation
                                                (heterogeneous P1), solid thermophysics and chemistry, porosity
tutorials/                                      test cases with reference results (see tutorials/README.md)
COPYING                                         GNU General Public License v3
```

## Running

A case has two regions, `fluidZone` (plenum, bed and freeboard, as in porousGasificationFoam) and `solidZone` (refractory). The region names are fixed. Run in serial with

```sh
porousGasificationFoamMR
```

or in parallel with `decomposePar -allRegions` and `mpirun -np N porousGasificationFoamMR -parallel`.

Note: in a multi-region case, the gas chemistry model reads `reactions` and `odeCoeffs` from the case-level `constant/chemistryProperties`, not from `constant/fluidZone/`. Changes written only in the region directory are silently ignored.

## Options added by this version

Unless stated otherwise, a feature is off, or reproduces the behaviour of porousGasificationFoam, when its keyword is absent.

| Feature | Where | Keyword (default) |
|---|---|---|
| Laminar EDM combustion | `constant/fluidZone/combustionProperties` | `combustionModel laminarEDM<psiReactionThermo>;` coefficients in `laminarEDMCoeffs`: `A` (4), `B` (0.5), `oxidant`, `reactions` (required) |
| Porosity-weighted realizability cap of laminarEDM | `laminarEDMCoeffs` | `porosityWeightedCap` (no) |
| Ergun resistance from the local porosity | cellZone `biomass_bed`; `system/controlDict` | always on if the zone exists; `ergunDp` particle diameter (0.05 m) |
| Inertial Ergun term | `system/controlDict` | `ergunInertialEps` (no); then set `forchheimerCoeff 0` in `constant/fluidZone/porosityProperties` |
| Kinetic energy update after the pressure correction | `system/controlDict` | `updateK` (yes since 1.1; no in 1.0) |
| Bed collapse | `constant/fluidZone/pyrolysisProperties`, `pyrolysisCoeffs` | `bedCollapse` (off), `shrinkA` (0.8), `collapseInterval` (20 s), `bedZoneName`, `phiS0Ref`, `phiS0Init`, `collapseMinPhis` (1e-3) |
| Held igniter | `pyrolysisCoeffs` | `igniterCellSet`, `igniterTemperature`, `igniterDuration` |
| Porosity ceiling from the ash | `pyrolysisCoeffs` | `ashDensity` (650 kg/m³) |
| Radiation–solid temperature sub-iterations | `constant/fluidZone/radiationProperties` | `nRadCorrectors` (1, off) |
| Under-relaxation of the wall radiative flux | `radiationProperties` | `qrRelax` (1, off; 0.3 used in the paper) |
| Time-step limits | `system/controlDict` | `maxRa`, `maxGa`, `maxUa` (0.3 each), `UaSmall` velocity floor of the `maxUa` criterion (1e-3 m/s); refractory diffusion number always on |
| Gas heat-release time-step limit | `constant/fluidZone/chemistryProperties` | `gasChemistryTimeStepControl` (false), `gasAlphaTemp` (0.05) |
| Bounded Arrhenius rate (temperature floored above 0 K) | `constant/reactions` | reaction type `irreversibleBoundedArrhenius` |

The defaults and the dictionaries above were checked against the `lookupOrDefault` calls of the source (22/10/2026).

### Corrections of version 1.1

Each correction of version 1.1 has a switch; with all of them set to `no`, version 1.1 reproduces version 1.0 bit for bit
(same log and fields on a 2 s restart of the reference case).

| Correction | Where | Keyword (default in 1.1) |
|---|---|---|
| Storage terms in conservative form, `ddt(porosityF, rho, X)` with the old-time porosity and density (species, energy), and the missing `p d(porosityF psi)/dt` part of the pressure equation; closes the global mass balance | `system/controlDict` | `conservativeDdt` (yes) |
| Boundary values of the bed indicator `whereIs` in the P1 model copied from the adjacent cell (they stayed at 1, so the wall diffusion coefficient included the bed absorption next to pure gas) | `radiationProperties` | `boundaryWhereIs` (yes) |
| Exact border factor of the bed: area of the faces shared with gas cells divided by the cell volume (the cell length used before gave 29% less on the bed top of a wedge mesh); written as `surfFrad` | `radiationProperties` | `exactBorderFactor` (yes) |
| Initial gas-chemistry time step applied at a start from t = 0 only, not at every restart, so that a restarted run continues the same sequence of time steps (write times with `writePrecision 17`) | `system/controlDict` | `exactRestart` (yes) |

## Test cases

See `tutorials/README.md`: the reference case of the paper, its variant with a lid, and a two-second restart with reference output and a comparison script (`kimReactor149_restart/check.sh`), which is the quick check of an installation.

## Citing

If you use this solver, please cite the accompanying article (S.K. Tchimoe et al., porousGasificationFoamMR: an OpenFOAM solver for transient fixed-bed biomass conversion with conjugate heat transfer and bed shrinkage, submitted to Computer Physics Communications, 2026) and the original porousGasificationFoam paper: P.J. Żuk, B. Tużnik, T. Rymarz, K. Kwiatkowski, M. Dudyński, F.C.C. Galeazzo, G.C. Krieger Filho, OpenFOAM solver for thermal and chemical conversion in porous media, Comput. Phys. Commun. 278 (2022) 108407, https://doi.org/10.1016/j.cpc.2022.108407.

## Licence

GNU General Public License v3 or later (see `COPYING`), as OpenFOAM and porousGasificationFoam.
