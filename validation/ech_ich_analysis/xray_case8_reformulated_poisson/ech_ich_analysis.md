# PICOS ECH/ICH Analysis

This analysis reads PICOS HDF5 output directly and uses the species charge sign to label ECH and ICH channels.

Moment convention:

- PICOS mesh moments are direct particle-weight sums over the sampled guiding-center ensemble.
- `P_perp = m n <v_perp^2>/2`; density and pressure deposition are not divided by `v_perp`.
- The `1/v_perp` factor appears only when reconstructing the full gyrotropic distribution from a `(v_parallel,v_perp)` histogram: `f = g/(2*pi*v_perp)`.

Summary rows:

- `xray_case8_reformulated_poisson` `spp_1` ICH: RF=0 Prf=300000.0 W, Tpar_mean=8.70181 eV, Tper_mean=20.7485 eV
- `xray_case8_reformulated_poisson` `spp_2` ECH: RF=0 Prf=300000.0 W, Tpar_mean=16.6757 eV, Tper_mean=17.3084 eV

Artifacts:

- `ech_ich_summary.csv`
- `ech_ich_profiles.png`
- `ech_ich_time_traces.png`
- `force_balance_spp_1.png`
- `velocity_space_spp_1_step1.png`
- `force_balance_spp_2.png`
- `velocity_space_spp_2_step1.png`
