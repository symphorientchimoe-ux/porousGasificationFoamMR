# porousGasificationFoamMR — user manual

This manual covers what differs from porousGasificationFoam. For the description of the porous medium, of the solid species and of the heterogeneous reactions (files `pyrolysisProperties`, `solidThermophysicalProperties`, `reactions`, `heatTransferProperties`, porosity fields), read the manual of porousGasificationFoam (`porousGasificationFoam_manual.pdf` in https://github.com/btuznik/porousGasificationFoam): those files keep their format. The equations are given in the accompanying paper.

## 1. Installation

```sh
source /path/to/OpenFOAM-8/etc/bashrc
./Allwmake
```

The six libraries go to `$FOAM_USER_LIBBIN`, the solver to `$FOAM_USER_APPBIN`. They have the same names as those of porousGasificationFoam: do not install both in the same user directory (see `README.md`).

Check the installation with the quick test (one core, a few minutes):

```sh
cd tutorials/kimReactor149_restart
./Allrun
./check.sh        # must end with "ALL CHECKS PASSED"
```

## 2. Structure of a case

A case has exactly two mesh regions with fixed names:

| Region | Content | Equations |
|---|---|---|
| `fluidZone` | plenum, fuel bed and freeboard, described as in porousGasificationFoam (porosity field, cellZone `biomass_bed`) | gas: mass, momentum, energy, species; solid inside the bed: energy, species, porosity |
| `solidZone` | refractory wall | heat conduction (standard OpenFOAM `heSolidThermo`) |

The two regions are usually obtained from one mesh with `splitMeshRegions -cellZones`. Each region has its own `constant/<region>/polyMesh`, `0/<region>` fields and `system/<region>/fvSchemes`, `fvSolution`. `constant/regionProperties` is not read: the names are fixed in the solver.

Files of the gas region (`constant/fluidZone/`): `thermophysicalProperties`, `thermo.compressibleGas`, `combustionProperties`, `chemistryProperties`, `radiationProperties`, `pyrolysisProperties`, `solidThermophysicalProperties`, `heatTransferProperties`, `porosityProperties`, `reactions`, `g`, `fvOptions` (optional). Files of the refractory (`constant/solidZone/`): `thermophysicalProperties` (heSolidThermo, constant properties).

**Pitfall.** The gas chemistry model reads `reactions` and `odeCoeffs` from the case-level `constant/chemistryProperties`, not from `constant/fluidZone/`. A change written only in the region file is silently ignored. The test cases avoid the problem with a case-level file that contains only `#include "fluidZone/chemistryProperties"`: keep that layout.

## 3. Coupling the two regions

On the interface patch of each region, the temperature uses the conjugate condition with radiation:

```
// 0/fluidZone/T
fluidZone_to_solidZone
{
    type            compressible::turbulentTemperatureRadCoupledMixed;
    Tnbr            T;
    qr              qr;
    qrNbr           none;
    value           $internalField;
}
// 0/solidZone/T
solidZone_to_fluidZone
{
    type            compressible::turbulentTemperatureRadCoupledMixed;
    Tnbr            T;
    qr              none;
    qrNbr           qr;
    value           $internalField;
}
// 0/fluidZone/G
fluidZone_to_solidZone
{
    type            MarshakRadiation;
    emissivityMode  lookup;
    emissivity      uniform 0.67;
    value           $internalField;
}
```

`qr` is the wall radiative flux computed by the P1 model on the gas side. The incident radiation `G` uses `MarshakRadiation` with the wall emissivity on the same patch. The outer face of the refractory can use `externalWallHeatFluxTemperature` (convection and radiation to the surroundings, optionally with a thin steel shell given by `thicknessLayers` and `kappaLayers`). Copy the complete conditions from `tutorials/kimReactor149_lid/60/`.

**Under-relaxation of qr.** The face temperature of the coupled wall depends on `qr`, and `qr` depends on that face temperature through the Marshak condition, one time step behind. On a hot freeboard wall the gain of this loop is close to one, and a bed collapse can push it past one: the face temperature then alternates from step to step until it diverges. Set `qrRelax 0.3;` in `constant/fluidZone/radiationProperties` (default 1, no relaxation). `qr` is written with the results and read back on restart, so a restarted run continues the same relaxation.

## 4. Gas-phase combustion: laminarEDM

```
// constant/fluidZone/combustionProperties
combustionModel  laminarEDM<psiReactionThermo>;
laminarEDMCoeffs
{
    A          4;        // Magnussen-Hjertager constant (default 4)
    B          0.5;      // product term constant (default 0.5)
    oxidant    O2;
    reactions             // fuel, oxidant-to-fuel mass ratio s, products
    (
        COoxidation  { fuel CO;  s 0.5714; products (CO2); }
        H2oxidation  { fuel H2;  s 7.937;  products (H2O); }
    );
}
```

The rate of each reaction is min(R_kin, A ρ|S| min(Y_f, Y_ox/s, B Y_p/(1+s))), with |S| = √2 |symm(∇U)| the resolved strain rate. The ratio κ = R/R_kin is written as the field `laminarEDM:kappa_<fuel>` (1: kinetics-limited, below 1: mixing-limited). The reaction names must be those of the gas `reactions` file; the tutorial `kimReactor149_lid` has four (CO, H2, CH4 and tar oxidation).

## 5. Bed resistance and grate

Inside the cellZone `biomass_bed`, the Darcy coefficient is recomputed at every time step from the local porosity, Df = 150(1−ε)²/(ε³ d_p²), with ε clamped to [0.05, 0.99]. The particle diameter is `ergunDp` in `system/controlDict` (default 0.05 m). With `ergunInertialEps yes;` the inertial Ergun term is added with the local porosity; then set `forchheimerCoeff 0` in `porosityProperties`.

A perforated grate can be represented by an `explicitPorositySource` (DarcyForchheimer, d = 0, F = ζ/L) in `constant/fluidZone/fvOptions`. **Spread it over at least five cell layers:** a single very resistive layer is bypassed by the face-flux interpolation and gives a pressure drop more than ten times too small (paper, Section 3.5). Always check the computed Δp against the target.

## 6. Bed collapse

```
// constant/fluidZone/pyrolysisProperties, pyrolysisCoeffs
bedCollapse       on;      // default off
shrinkA           0.8;     // Yang's shrinkage coefficient a
collapseInterval  20;      // s, time between two collapses
phiS0Ref          0.62;    // reference solid fraction (default: from the initial field)
bedZoneName       biomass_bed;
```

Every `collapseInterval`, each vertical column of the bed is shortened according to Yang's law, h = Δz (φ_s0 − a(φ_s0 − φ_s))/φ_s0Ref, and the solid is moved down cell by cell with exact conservation of the mass of each solid species. The bed must be a structured column/layer mesh with a vertical axis. The log reports each collapse:

```
bedCollapse at t = 2400.5 (shrinkA 0.8, phiS0Ref 0.62)
    bed top axis column : 0.35 -> 0.345001255263 m ; all columns min/max : ...
    cells emptied 1, top residues merged 1, h clamped to dz 6, porosity<0 clipped 0, ...
    conservation M(Yfuel) [kg] : before ... after ... rel 2.39e-16
```

## 7. Igniter

```
// pyrolysisCoeffs
igniterCellSet       igniterZone;   // a cellSet in constant/fluidZone/polyMesh/sets
igniterTemperature   1000;          // K
igniterDuration      300;           // s
```

The solid temperature of the igniter cells is held at `igniterTemperature` until `igniterDuration`. The values are read only at start-up: to change them during a run, stop it, edit the file and restart (see `tutorials/kimReactor149/Allrun`).

## 8. Time step

On top of the Courant (`maxCo`) and solid diffusion (`maxDi`) numbers of porousGasificationFoam, the time step is bounded by the refractory diffusion number (always on), and by three relaxation numbers in `system/controlDict`: `maxRa` (radiative sink), `maxGa` (gas–solid exchange) and `maxUa` (velocity growth), 0.3 by default. With `gasChemistryTimeStepControl true;` and `gasAlphaTemp` (0.05) in `constant/fluidZone/chemistryProperties`, the gas heat release limits the relative temperature change per step. Each step the log line `dtLimits Co ... Di ... Ra ... Ga ... Ua ...` gives, for each criterion, the ratio of its largest allowed step to the current step: the smallest ratio is the limiting criterion.

## 9. Linear solvers

Use `PBiCGStab` (with DILU) for the non-symmetric equations of the gas region (U, h, species, solid species, porosity). With `PBiCG`, the runs of the paper eventually failed after a calm period, with a sudden jump of the gas temperature (paper, Section 3.6). Keep the residual output enabled: a failed linear solve is otherwise invisible.

## 10. Outputs and diagnostics

Fields added to the usual ones: `qr` (wall radiative flux, W m⁻²), `laminarEDM:kappa_<fuel>`, and the solid fields of porousGasificationFoam (`Ts`, `porosityF`, solid mass fractions). Each step the solver checks the main fields and matrices for non-finite values; at the first one it stops and prints the stage and the first cells concerned with their context (lines starting with `DBG`). Finite values outside the physical range are only reported.

## 11. Running long cases

Run in parallel with `decomposePar -allRegions` and `mpirun -np N porousGasificationFoamMR -parallel`. When several runs share a machine, pin each one to its own cores (`taskset -c ... mpirun --bind-to none ...`): by default mpirun binds every run to the same cores.

A restart from a written time is exact only if every lagged quantity is written. Use `writePrecision` 12 or more (or binary format). Restarts of the runs of the paper changed the temperatures above the front by up to about 200 K later in the run (paper, Section 5.5): compare runs restarted from the same state, not a restarted run with an uninterrupted one.
