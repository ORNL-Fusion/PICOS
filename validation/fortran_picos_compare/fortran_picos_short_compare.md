# Fortran versus PICOS++ Short Comparison

Fortran case: `/Users/78k/Desktop/picos_kinetic_electron_eval/LinearFokkerPlanck_Axisymmetric/InputFiles/xp_Case8PicosCompare.in`
PICOS++ tags: `xray_case8_fortran_compare_nonrel`, `xray_case8_fortran_compare_rel`

| run | N | mean E [eV] | P95 [eV] | P99 [eV] | max E [eV] | mean z [m] | mean pitch | max v/c |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| fortran_nonrel | 640 | 23.93 | 70.22 | 125.5 | 252.5 | -0.01799 | -0.02093 | 0.03144 |
| picos_nonrel | 640 | 28.42 | 64.08 | 83.2 | 101.4 | -1.569 | -0.02734 | 0.01992 |
| picos_rel | 640 | 28.42 | 66.12 | 83.44 | 101.4 | -1.538 | -0.02043 | 0.01992 |

PICOS++ nonrel mean-energy delta from Fortran: `4.48704 eV`.
PICOS++ relativistic mean-energy delta from Fortran: `4.48878 eV`.

Validation plots:

- `ech_operator_energy_hist.png`
- `ech_operator_energy_cdf.png`
- `ech_operator_energy_vs_z.png`
- `ech_operator_energy_vs_pitch.png`
- `ech_operator_velocity_space.png`
- `ech_operator_energy_metrics.png`

Interpretation:

- PICOS++ nonrel mean energy is `4.48704 eV` higher than Fortran.
- PICOS++ nonrel P99 energy is `-42.3269 eV` different from Fortran.
- Fortran produces a hotter extreme tail here: max energy `252.518 eV` versus PICOS++ nonrel `101.353 eV`.
- The comparison is stochastic and not particle-by-particle matched; it validates ensemble behavior of the ECH operator.
