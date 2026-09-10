# PICOS++ HDF5 Restart Workflow

This restart path is intended for two-stage kinetic ion/electron MPEX runs:

1. Run plasma formation to a quasi-steady state with source, collisions, and field solve enabled, but RF/ECH off.
2. Start a fresh output tag, load the last HDF5 snapshot from stage 1, then enable RF/ECH.

The restart loader runs after normal input allocation and before characteristic normalization. It reads SI-unit HDF5 output data and then lets the normal PICOS++ normalization and PIC constructor rebuild particle cells, particle fields, and deposited moments.

Use the same `mpisForFields` and total MPI rank count in stage 2 as in stage 1. The restart reader loads rank-local `FIELDS_FILE_<rank>.h5` and `PARTICLES_FILE_<rank>.h5` files; it does not repartition particles between a different MPI decomposition.

## Required Stage-1 Output

Set every species to write all particles:

```text
pctSupPartOutput1             100
pctSupPartOutput2             100
```

The restart requires these particle datasets in every `PARTICLES_FILE_<rank>.h5`:

```text
X_p
V_p
a_p
```

`mu_p` is restored when present. Mesh fields are restored when the corresponding field datasets are present; otherwise the run continues with fields initialized from the input profiles and the enabled field solver updates them.

## Input Keys

Use these keys in the stage-2 input file:

```text
restart_enabled             1
restart_path                /path/to/stage1/output/HDF5
restart_snapshot            -1
restart_continueTime        0
restart_particleFilePrefix  PARTICLES_FILE_
restart_fieldsFilePrefix    FIELDS_FILE_
```

`restart_snapshot -1` loads the latest numeric snapshot. Use a non-negative integer to load a specific snapshot.

`restart_continueTime 0` resets the stage-2 clock to zero after loading the particle state. This is usually best for ECH turn-on because `RF_*_t_ON` and `RF_*_t_OFF` then refer to the stage-2 time window.

`restart_continueTime 1` continues from the saved physical time in the restart file. Use this only when RF/source timing in the stage-2 input is written on the absolute clock of the combined run.

## MPEX Stage Pattern

For the scenario-14 MPEX sequence currently used here, the steady-state tag is:

```text
mpex_scenario14_ex8_steady_2ms_p262144_coll_mpexprof_nersc_nonrel
```

It runs for `2.0e-3` s with `40` output intervals, so the saved cadence is `50` microseconds plus the initial snapshot.

Stage 1, steady state:

```text
SW_EfieldSolve              1
SW_fieldSolveModel          2
SW_Collisions               1
SW_pairSource               1
SW_RFheating                0
restart_enabled             0
```

Run with a steady-state tag:

```bash
mpirun -np 128 ${PICOS_BIN} 1-D outputFiles mpex_scenario14_ex8_steady_2ms_p262144_coll_mpexprof_nersc_nonrel
```

Stage 2, ECH continuation:

```text
SW_EfieldSolve              1
SW_fieldSolveModel          2
SW_Collisions               1
SW_pairSource               1
SW_RFheating                1
SW_RFheatingIons            0
SW_RFheatingElectrons       1
restart_enabled             1
restart_path                outputFiles/mpex_scenario14_ex8_steady_2ms_p262144_coll_mpexprof_nersc_nonrel/HDF5
restart_snapshot            -1
restart_continueTime        0
```

Run with a new output tag:

```bash
mpirun -np 128 ${PICOS_BIN} 1-D outputFiles mpex_s14_ech_from_steady_nonrel
```

Do not restart into the same output tag as stage 1, because the HDF5 writer truncates snapshot files at the beginning of a new run.

## Scenario-14 Deck Generation

`scripts/compare_mpex_scenario14_fig17.py` supports the restart keys directly:

```bash
python3 scripts/compare_mpex_scenario14_fig17.py \
  --setup \
  --picos-run-mode nonrel \
  --picos-tag-nonrel mpex_s14_well_steady_nonrel \
  --rf-heating 0 \
  --collisions 1 \
  --efield-solve 1 \
  --field-solve-model 2 \
  --reformulated-poisson-quasineutral 1 \
  --pair-source 1

python3 scripts/compare_mpex_scenario14_fig17.py \
  --setup \
  --picos-run-mode nonrel \
  --picos-tag-nonrel mpex_s14_well_ech_restart_nonrel \
  --rf-heating 1 \
  --restart-path /path/to/outputFiles/mpex_s14_well_steady_nonrel/HDF5 \
  --restart-snapshot -1 \
  --collisions 1 \
  --efield-solve 1 \
  --field-solve-model 2 \
  --reformulated-poisson-quasineutral 1 \
  --pair-source 1
```
