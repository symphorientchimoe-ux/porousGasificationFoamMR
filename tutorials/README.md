# Test cases

| Case | What it runs | Cores | Run time |
|---|---|---|---|
| `kimReactor149_lid` | Reference validation run of the paper (Table 7): counter-current wood-pellet bed of Kim et al. (2016), primary air 149 kg m⁻² h⁻¹, 5 mm bed cells, reactor closed by a lid; starts from the open-top state at 60 s | 2 | about 5 days |
| `kimReactor149` | Same reactor with an open top, from ignition (t = 0): the comparison run of the paper, and the source of the 60 s state used by `kimReactor149_lid` | 2 | about 35 h |
| `kimReactor149_restart` | Two seconds of the reference run (`kimReactor149_lid`), restarted from its state at 2400 s, with a bed collapse every 0.5 s. Every model is active: pyrolysis and char reactions, gas combustion (laminarEDM), P1 radiation, refractory conduction, bed collapse. Reference output included | 1 | a few minutes |

`kimReactor149_restart` is the quick check of an installation. Run `./Allrun`, then `./check.sh`, which compares the probe values, the bed heights and the conservation lines of the log with the reference files in `reference/`.

## Origin of the input data

The geometry, the fuel and the operating conditions come from Kim et al., Waste Management 49 (2016). The physical parameters, with their sources, are listed in the accompanying paper (Table 2 for the physical parameters, Table 3 for the numerical settings). The igniter temperature and duration are adjusted parameters (the paper describes how).
