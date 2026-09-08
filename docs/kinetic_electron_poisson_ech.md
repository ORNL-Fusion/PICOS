# Kinetic-Electron Poisson/ECH Prototype

This branch keeps the original 1D-2V guiding-center PICOS workflow as the recommended baseline and adds a selectable electrostatic Poisson field solve, kinetic negative-charge species, and an optional 1D-3V Boris pusher for development tests.

## Field Solver Modes

Use these input flags in `inputFiles/input_file*.input`:

```text
SW_EfieldSolve              1
SW_fieldSolveModel          0   // original quasi-neutral Ohm-law hybrid solve
SW_fieldSolveModel          1   // electrostatic Poisson solve from kinetic charge density
SW_fieldSolveModel          2   // reformulated Poisson electric-field update from kinetic stress moments
```

`SW_fieldSolveModel` is optional. If it is absent, PICOS uses the original Ohm-law model.

For Poisson runs, the optional voltage boundary values are:

```text
IC_phiLeft                  0.0
IC_phiRight                 0.0
```

They are interpreted as physical volts in the input file and normalized internally by `E0*L0`.

Poisson boundary models are:

```text
Poisson_BCModel             0   // Dirichlet phiLeft/phiRight
Poisson_BCModel             1   // periodic electrostatic potential
Poisson_BCModel             2   // simple sheath-shifted wall potential
Poisson_sheathCoefficient   3.0
```

For `Poisson_BCModel=2`, low parallel-energy kinetic electrons reflect from the sheath barrier instead of immediately being counted as boundary losses.

For the reformulated Poisson model from `/Users/78k/Downloads/Reformulated Poissons equation.pdf`, add:

```text
ReformulatedPoisson_lambda       -1.0   // physical meters if positive; -1 uses the code's normalized Poisson scale
ReformulatedPoisson_quasiNeutral 0      // 1 uses E = div(S_i-S_e)/(n_i+n_e/epsilon)
```

The implemented 1D form is:

```text
dE_x/dt + omega_p^2 E_x = d(S_i-S_e)/dx / lambda^2
omega_p^2 = (n_i + n_e/epsilon) / lambda^2
epsilon = m_e / m_i
```

`S_i-S_e` is computed from the kinetic particle moments on the mesh. PICOS stores `P11_m` as a mass-weighted parallel second moment, so the solver uses `P11_m/M` for each self-consistent positive/negative species before taking the 1D divergence. The updated electric field is advanced implicitly in the local plasma-frequency term, and `Phi_m` is reconstructed from `E_x = -d(phi)/dx` for output and sheath diagnostics.

The ready-to-edit template is `templateFILES/input_file_reformulated_poisson.input`.

## Kinetic Electrons

Add electrons as a self-consistent particle species in `ions_properties*.ion`:

```text
SPECIES2                    1
Z2                          -1
M2                          0.000548579909
```

`M2` is the electron mass in atomic mass units. The code still uses the historical `IONS` container name, but the species handling now permits negative-charge kinetic species.

Use `advanceParticleMethod 1` for the guiding-center baseline. In this mode ions and electrons are kinetic particles, but `V_p` keeps the historical two-column `v_parallel, v_perp` storage. Use `advanceParticleMethod 3` only when explicitly testing the optional 1D-3V Boris/full-orbit prototype; in that mode `V_p` is written as three velocity components.

## Poisson Mapping

After particles deposit moments, the Poisson solver assembles normalized charge density on the field ranks:

```text
rho = sum_s (q_s/q0) * n_s / (n0 * V0)
```

Only `numberOfParticleSpecies` self-consistent species are included. Tracer species are excluded from the field charge density.

The solver then computes:

```text
d2(phi)/dx2 = -rho / epsilon0
E_x = -d(phi)/dx
```

The solved potential can be written to HDF5 by adding `Phi_m` to `outputs_variables`. `Phi_m` is also exchanged with particle ranks for sheath boundary handling.

The particle moment deposition now folds boundary ghost-cell shape-function contributions back into the physical mesh. Periodic species fold left/right support to the opposite side; nonperiodic species fold support into the nearest wall cell. This keeps deposited density/charge conserved for particles whose shape function straddles a boundary.

## Moment Deposition Convention

For the recommended `advanceParticleMethod 1` guiding-center model, `V_p` stores `(v_parallel, v_perp)`, where `v_perp` is the perpendicular speed magnitude. Particles are sampled from the reduced 1D-2V guiding-center distribution, so the cylindrical velocity-space Jacobian is already represented by the particle ensemble. Mesh moments are therefore direct particle-weight sums:

```text
n       = sum_p w_p
n u     = sum_p w_p v_parallel
P_11    = sum_p w_p m v_parallel^2
P_perp  = sum_p w_p m v_perp^2 / 2
```

Do not divide mesh density, flux, or pressure deposition by `v_perp`. A `1/v_perp` factor is only relevant when reconstructing a full distribution function from a histogram in `(v_parallel, v_perp)`; in that diagnostic case the bin normalization must include the cylindrical measure `2*pi*v_perp*dv_parallel*dv_perp`.

## Electron Time-Step Resolution

When a kinetic electron species is present, the time-step selection now considers:

```text
all-particle CFL:             DX / max(|v|)
electron CFL:                 DX / max(|v_e|)
electron gyro period:         DTc * (2*pi/Omega_ce)
electron plasma time scale:   DTc * (1/omega_pe)
```

The electron plasma limiter uses `1/omega_pe`, not a full oscillation period.

## Relativistic Electron Correction

For high-energy ECH tests, enable:

```text
SW_relativisticElectrons    1
```

This currently applies to negative-charge kinetic species only. It does three things:

```text
1. Uses gamma-corrected cyclotron resonance: Omega_c -> Omega_c/gamma.
2. Applies RF/ECH energy kicks using KE = (gamma - 1) m c^2.
3. Caps post-push electron speed below c for relativistic-electron runs.
```

This is still a 1D-2V guiding-center model. It is not yet a full relativistic particle pusher with momentum as the evolved variable, so strong electrostatic acceleration and sheath validation still need dedicated tests.

## ECH/RF Notes

The RF operator supports negative-charge resonances. For electron cyclotron heating use approximately:

```text
SW_RFheating                1
SW_RFheatingIons            0
SW_RFheatingElectrons       1
SW_relativisticElectrons    1
RF_electron_n_harmonic      1
RF_electron_freq            e*B/(2*pi*m_e)
RF_electron_handedness      -1
```

`SW_RFheating` is the master RF switch. `SW_RFheatingIons` gates RF kicks for positive-charge species, and `SW_RFheatingElectrons` gates RF kicks for negative-charge species. Set both to `1` for a mixed ion/electron RF test, or set only one to isolate the heating channel.

RF/ECH settings are species-specific. Use `RF_ion_*` for positive-charge species and `RF_electron_*` for negative-charge species:

```text
RF_ion_Prf                      50E3
RF_ion_n_harmonic               1
RF_ion_freq                     8.385E6
RF_ion_x1                       4.0
RF_ion_x2                       6.5
RF_ion_t_ON                     12000
RF_ion_t_OFF                    20000
RF_ion_kpar                     20
RF_ion_kper                     100
RF_ion_handedness               -1
RF_ion_EfieldMode               0
RF_ion_EfieldAmplitude          0.0
RF_ion_maxEnergyGainFraction    0.0
RF_ion_maxParticleEnergy        0.0
RF_ion_maxVelocityFractionC     0.0
RF_ion_Prf_fileName             Prf_profile.txt
RF_ion_Prf_NS                   200

RF_electron_Prf                      3.0E5
RF_electron_n_harmonic               2
RF_electron_freq                     6.5E10
RF_electron_x1                       -0.4
RF_electron_x2                       0.0
RF_electron_t_ON                     0.0
RF_electron_t_OFF                    2.0E-3
RF_electron_kpar                     5.864E3
RF_electron_kper                     2.992E4
RF_electron_handedness               -1
RF_electron_EfieldMode               0
RF_electron_EfieldAmplitude          1.0E4
RF_electron_maxEnergyGainFraction    0.25
RF_electron_maxParticleEnergy        5000
RF_electron_maxVelocityFractionC     0.2
RF_electron_Prf_fileName             Prf_profile.txt
RF_electron_Prf_NS                   200
```

Legacy `RF_*` keys are still accepted. When present, they are used as defaults for both species unless a species-specific key overrides them.

The operator now keeps electron mass as floating point and skips RF kicks safely when no particles are resonant during a time step. It also supports fixed archived RF electric field amplitude:

```text
RF_electron_EfieldMode       0   // absorbed-power balance using RF_electron_Prf
RF_electron_EfieldMode       1   // fixed RF_electron_EfieldAmplitude [V/m]
RF_electron_EfieldAmplitude  10000
```

For kinetic-electron prototype runs, use the optional nonrelativistic guard rails:

```text
RF_electron_maxEnergyGainFraction    0.25
RF_electron_maxParticleEnergy        5000      // eV; 0 disables
RF_electron_maxVelocityFractionC     0.2       // fraction of c; 0 disables
```

These are especially useful with archived `Ew` inputs before a fully validated ECH operator is available.

## X-ray/ECH Benchmark

The archived 2020 Proto-MPEX X-ray study lives at:

```text
/Users/78k/Desktop/2020_07_13_XrayStudy
```

Generate the archived-data summary and plots:

```bash
/usr/bin/env PYTHONPYCACHEPREFIX=/private/tmp/picos_pycache MPLCONFIGDIR=/private/tmp/picos_mpl \
  python3 scripts/compare_ech_xray_study.py --cases Case8 Case11 Case15 --out-dir validation/xray_study
```

Outputs:

```text
validation/xray_study/ech_xray_case_summary.csv
validation/xray_study/ech_xray_summary.md
validation/xray_study/ech_xray_final_energy_hist.png
validation/xray_study/ech_xray_energy_vs_z.png
```

## ECH/ICH HDF5 Analysis

Use the repo-local HDF5 analysis scripts instead of editing hard-coded `myDir` paths in the older ProtoLite MATLAB/Python files:

```bash
PYTHONPYCACHEPREFIX=/private/tmp/picos_pycache MPLCONFIGDIR=/private/tmp/picos_mpl \
  python scripts/analyze_picos_ech_ich.py xray_case8_reformulated_poisson \
  --out-dir validation/ech_ich_analysis/xray_case8_reformulated_poisson \
  --velocity-hist
```

The script accepts a run tag, a run directory, or a direct `HDF5` directory. It reads `/rf/ion` and `/rf/electron` metadata from `main.h5`, labels positive-`Z` species as ICH-channel species and negative-`Z` species as ECH-channel species, and writes:

```text
ech_ich_summary.csv
ech_ich_analysis.md
ech_ich_profiles.png
ech_ich_time_traces.png
force_balance_spp_*.png
velocity_space_spp_*_step*.png
```

The velocity-space diagnostic writes the reduced guiding-center distribution and the full gyrotropic distribution using `f = g/(2*pi*v_perp)`. This is only a diagnostic reconstruction; it is not used in PIC mesh moment deposition.

Current archived-data reference:

```text
Case8 heated:     mean final E = 29.7 eV, P95 = 65.9 eV, P99 = 161 eV, >10Te fraction = 0.0101
Case11 no heat:   mean final E = 24.3 eV, P95 = 62.4 eV, P99 = 90.7 eV, >10Te fraction = 0
Bres for 65 GHz, harmonic 2: 1.161 T
```

Create PICOS decks from the archived cases:

```bash
python3 scripts/create_picos_xray_case.py --case Case8 --tag xray_case8_gc_direct_fullsetup --field-solve none --advance-particle-method 1
python3 scripts/create_picos_xray_case.py --case Case11 --tag xray_case11_gc_direct_fullsetup --field-solve none --advance-particle-method 1
```

Use `--field-solve none` for direct comparison to the archived X-ray cases because those namelists use `iPotential = .false.`. High-density Case8/Case11 are not good Poisson smoke tests unless the mesh is Debye-resolved; the current coarse smoke had `lambda_D/DX = 2.5e-5` and Poisson noise drove unphysical velocities.

Cleaned Case8 input/template deck set:

```text
input_file_xray_case8_fortran_compare_nonrel.input     full Case8 ECH comparison, field solve off, nonrelativistic electrons
input_file_xray_case8_fortran_compare_rel.input        full Case8 ECH comparison, field solve off, relativistic electrons
input_file_xray_case8_gc_rf_smoke.input                short RF/ECH smoke, species-split RF_electron_* keys, field solve off
input_file_xray_case8_reformulated_poisson.input       short reformulated-Poisson smoke, RF off, SW_fieldSolveModel=2
```

The same named `ions_properties_*.ion` files and profile files are kept under `picosFILES/inputFiles` for local runs and mirrored under `templateFILES` for the branch.

Short direct RF smoke deck:

```bash
python3 scripts/create_picos_xray_case.py --case Case8 --tag xray_case8_gc_rf_smoke \
  --physical-time 1.0e-12 --nx 80 --ion-particles-per-cell 2 --electron-particles-per-cell 8 \
  --output-saves 2 --collisions 0 --field-solve none --advance-particle-method 1 \
  --rf-heat-ions 0 --rf-heat-electrons 1 --relativistic-electrons 1

cd picosFILES
HDF5_USE_FILE_LOCKING=FALSE OMP_NUM_THREADS=1 mpirun -np 4 ../build/picosFILES/src/xpicos 1-D "$PWD/outputFiles" xray_case8_gc_rf_smoke
```

Run the reformulated-Poisson smoke deck:

```bash
cd picosFILES
HDF5_USE_FILE_LOCKING=FALSE OMP_NUM_THREADS=1 mpirun -np 4 ../build/picosFILES/src/xpicos 1-D "$PWD/outputFiles" xray_case8_reformulated_poisson
```

Summarize PICOS particle output without `h5py`:

```bash
python3 scripts/summarize_picos_particles_h5.py \
  picosFILES/outputFiles/xray_case8_gc_rf_smoke/HDF5/PARTICLES_FILE_0.h5 --iteration 1 --species spp_2 --relativistic
```

The guiding-center smoke output should report `V_p` with shape `{2,N}`. Current direct RF smoke result:

```text
finite = true
V_p shape = [2, 80]
superluminal particles = 0
max speed = 0.0181 c
mean relativistic energy = 26.5 eV
P95 relativistic energy = 68.0 eV
max relativistic energy = 83.8 eV
```

## Smoke Test

The branch includes:

```text
templateFILES/input_file_fully_kinetic.input
templateFILES/ions_properties_fully_kinetic.ion
```

To run from the repository root after copying templates into `picosFILES/inputFiles`:

```bash
mkdir -p picosFILES/inputFiles picosFILES/outputFiles
cp templateFILES/input_file_fully_kinetic.input picosFILES/inputFiles/
cp templateFILES/ions_properties_fully_kinetic.ion picosFILES/inputFiles/
cp templateFILES/MPEX_B_norm_PICOS_scenario_14.txt picosFILES/inputFiles/
cp templateFILES/MPEX_Tpar_norm_scenario_14.txt picosFILES/inputFiles/
cp templateFILES/MPEX_Tper_norm_scenario_14.txt picosFILES/inputFiles/
cp templateFILES/MPEX_n_norm_scenario_14.txt picosFILES/inputFiles/
cp templateFILES/ProtoMPEX_Te_norm_PICOS_c.txt picosFILES/inputFiles/
cd picosFILES
HDF5_USE_FILE_LOCKING=FALSE OMP_NUM_THREADS=1 mpirun -np 4 ../build/picosFILES/src/xpicos 1-D "$PWD/outputFiles" fully_kinetic
```

Expected HDF5 outputs include:

```text
HDF5/main.h5:/fieldSolveModel
HDF5/main.h5:/relativisticElectrons
HDF5/main.h5:/rf/heatIons
HDF5/main.h5:/rf/heatElectrons
HDF5/FIELDS_FILE_*.h5:/*/fields/EX_m/x
HDF5/FIELDS_FILE_*.h5:/*/fields/Phi_m/x
HDF5/PARTICLES_FILE_*.h5:/*/ions/spp_1
HDF5/PARTICLES_FILE_*.h5:/*/ions/spp_2
```

For this template, `advanceParticleMethod 1` is intentional. It is a kinetic-ion plus kinetic-electron 1D-2V guiding-center run, not the optional 1D-3V Boris prototype.

## Algorithm Checks

Run lightweight invariants:

```bash
python3 scripts/check_pic_algorithm_invariants.py --particles 10000 --boris-steps 2000
```

Current result:

```text
TSC weight-sum error:        2.22e-16
periodic deposited total:    conserved exactly for 10000 particles
wall deposited total:        conserved exactly for 10000 particles
Boris energy error:          2.51e-15 relative over 2000 no-E uniform-B steps
```

The Boris check covers the optional method-3 path. The production-style baseline remains method 1 unless a case explicitly asks for 1D-3V testing.

## Remaining Production Work

This is now a working electrostatic kinetic-electron prototype inside the existing 1D PICOS structure, not a complete production fully kinetic plasma code yet.

Main next items:

```text
1. Validate the ECH kick operator quantitatively against the 2020 distribution evolution, not only final-energy smoke statistics.
2. Add a Debye-resolution abort or explicit override for high-density Poisson runs.
3. Add charge-conserving current deposition if moving beyond electrostatic 1D.
4. Add real open-field Poisson/sheath boundary conditions tied to wall currents.
5. Keep the 1D-2V guiding-center path as the validated baseline while deciding whether any target physics really requires the optional 1D-3V or multidimensional full-orbit extension.
```
