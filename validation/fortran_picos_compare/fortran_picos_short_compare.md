# Fortran versus PICOS++ Short Comparison

Fortran case: `/Users/78k/Desktop/picos_kinetic_electron_eval/LinearFokkerPlanck_Axisymmetric/InputFiles/xp_Case8PicosCompare.in`
PICOS++ tags: `xray_case8_fortran_compare_nonrel`, `xray_case8_fortran_compare_rel`

| run | N | mean E [eV] | P95 [eV] | P99 [eV] | max E [eV] | mean z [m] | mean pitch | max v/c |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| fortran_nonrel | 640 | 24.27 | 77.3 | 128 | 241.6 | 0.06218 | -0.001925 | 0.03075 |
| picos_nonrel | 640 | 29.61 | 67.61 | 97.15 | 130.1 | -1.504 | 0.02544 | 0.02257 |
| picos_rel | 640 | 30.67 | 71.71 | 96.63 | 121.5 | -1.471 | 0.03047 | 0.02181 |

PICOS++ nonrel mean-energy delta from Fortran: `5.33635 eV`.
PICOS++ relativistic mean-energy delta from Fortran: `6.39387 eV`.
