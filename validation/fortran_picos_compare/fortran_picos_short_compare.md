# Fortran versus PICOS++ Short Comparison

Fortran case: `/Users/78k/Desktop/picos_kinetic_electron_eval/LinearFokkerPlanck_Axisymmetric/InputFiles/xp_Case8PicosCompare.in`
PICOS++ tags: `xray_case8_fortran_compare_nonrel`, `xray_case8_fortran_compare_rel`

| run | N | mean E [eV] | P95 [eV] | P99 [eV] | max E [eV] | mean z [m] | mean pitch | max v/c |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| fortran_nonrel | 640 | 24.27 | 77.3 | 128 | 241.6 | 0.06218 | -0.001925 | 0.03075 |
| picos_nonrel | 640 | 29.29 | 64.74 | 81.89 | 119.6 | -1.523 | -0.006839 | 0.02164 |
| picos_rel | 640 | 30.66 | 74.63 | 97.66 | 126.3 | -1.374 | -0.002256 | 0.02223 |

PICOS++ nonrel mean-energy delta from Fortran: `5.01533 eV`.
PICOS++ relativistic mean-energy delta from Fortran: `6.38825 eV`.

Validation plots:

- `ech_operator_energy_hist.png`
- `ech_operator_energy_cdf.png`
- `ech_operator_energy_vs_z.png`
- `ech_operator_energy_vs_pitch.png`
- `ech_operator_velocity_space.png`
- `ech_operator_energy_metrics.png`

Interpretation:

- PICOS++ nonrel mean energy is `5.01533 eV` higher than Fortran.
- PICOS++ nonrel P99 energy is `-46.0827 eV` different from Fortran.
- Fortran produces a hotter extreme tail here: max energy `241.598 eV` versus PICOS++ nonrel `119.612 eV`.
- The comparison is stochastic and not particle-by-particle matched; it validates ensemble behavior of the ECH operator.
