# Test cases

| Case | What it runs | Cores | Run time |
|---|---|---|---|
| `kimReactor149` | Reference validation run of the paper (version 1.1, Section 5, Table 7): counter-current wood-pellet bed of Kim et al. (2016), primary air 149 kg m⁻² h⁻¹, 5 mm bed cells, open top, from ignition (t = 0) | 2 | about 70 simulated s per hour; front at T1 after about 4700 s |
| `kimReactor149_lid` | Same reactor closed by a lid with a central outlet, as in the experiment; starts from the version 1.1 state of `kimReactor149` at 60 s. Its version 1.0 run is the "v1.0, lid" column of Table 7 | 2 | about three times slower than `kimReactor149` |
| `kimReactor149_restart` | Two seconds of the reference run, restarted from its version 1.1 state at 2688 s, with a bed collapse every 0.5 s. Every model is active: pyrolysis and char reactions, gas combustion (laminarEDM), P1 radiation, refractory conduction, bed collapse. Reference output included | 1 | a few minutes |

`kimReactor149_restart` is the quick check of an installation. Run `./Allrun`, then `./check.sh`, which compares the probe values, the bed heights and the conservation lines of the log with the reference files in `reference/`.

## Origin of the input data

The geometry, the fuel and the operating conditions come from Kim et al., Waste Management 49 (2016). The physical parameters, with their sources, are listed in the accompanying paper (Table 2 for the physical parameters, Table 3 for the numerical settings). The igniter temperature and duration are adjusted parameters (the paper describes how).
