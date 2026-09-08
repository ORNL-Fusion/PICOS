# Fortran versus PICOS++ Short Comparison

Fortran case: `/Users/78k/Desktop/picos_kinetic_electron_eval/LinearFokkerPlanck_Axisymmetric/InputFiles/xp_Case8PicosCompareTail.in`
PICOS++ tags: `xray_case8_fortran_compare_tail_nonrel`, `xray_case8_fortran_compare_tail_rel`

| run | N | mean E [eV] | P95 [eV] | P99 [eV] | max E [eV] | mean z [m] | mean pitch | max v/c |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| fortran_nonrel | 51200 | 24 | 74.15 | 117.9 | 273.9 | 0.004635 | -0.001872 | 0.03274 |
| picos_nonrel | 51200 | 24.12 | 74.55 | 117.4 | 308.2 | 0.02022 | -0.001023 | 0.03473 |
| picos_rel | 51200 | 24 | 74.5 | 119.3 | 335.5 | -0.01633 | -0.00196 | 0.03622 |

PICOS++ nonrel mean-energy delta from Fortran: `0.116353 eV`.
PICOS++ relativistic mean-energy delta from Fortran: `-0.00496224 eV`.

Validation plots:

- `ech_operator_energy_hist.png`
- `ech_operator_energy_cdf.png`
- `ech_operator_energy_vs_z.png`
- `ech_operator_energy_vs_pitch.png`
- `ech_operator_velocity_space.png`
- `ech_operator_energy_metrics.png`

Tail diagnostics:

| run | P99.5 [eV] | P99.9 [eV] | max E [eV] | frac E>5Te | frac E>10Te | frac E>15Te |
|---|---:|---:|---:|---:|---:|---:|
| fortran_nonrel | 136.3 | 189.8 | 273.9 | 0.0405664 | 0.00230469 | 9.76563e-05 |
| picos_nonrel | 136.2 | 191.8 | 308.2 | 0.0403125 | 0.00228516 | 0.000175781 |
| picos_rel | 141 | 185.1 | 335.5 | 0.0410352 | 0.00253906 | 0.000136719 |

Interpretation:

- PICOS++ nonrel mean energy differs from Fortran by `0.116353 eV`.
- PICOS++ nonrel P99 energy is `-0.472664 eV` different from Fortran.
- The single-particle maximum differs by `34.285 eV`; because this is stochastic, use P99/P99.5/P99.9 and threshold fractions for tail validation.
- The comparison is stochastic and not particle-by-particle matched; it validates ensemble behavior of the ECH operator.

RF diagnostics:

| run | RF event rate [1/s] | absorbed RF power [W] | Erf [V/m] | uE3 [W/(V/m)^2] |
|---|---:|---:|---:|---:|
| fortran_nonrel | 2.20948e+30 | -3.2668e+07 | nan | nan |
| picos_nonrel | nan | 2.79708e+08 | 0.882086 | 5.85075e+26 |
| picos_rel | nan | 1.88337e+08 | 0.280296 | 5.79428e+27 |

Validation checks:

| check | status | value | limit |
|---|---|---:|---:|
| nonrel particle count | PASS | 0 | 0 |
| rel particle count | PASS | 0 | 0 |
| nonrel mean-energy relative delta | PASS | 0.00484773 | 0.05 |
| rel mean-energy relative delta | PASS | 0.000206747 | 0.05 |
| nonrel P95-energy relative delta | PASS | 0.00543284 | 0.25 |
| rel P95-energy relative delta | PASS | 0.00472372 | 0.25 |
| nonrel P99-energy relative delta | PASS | 0.00400937 | 0.35 |
| rel P99-energy relative delta | PASS | 0.0117409 | 0.35 |
| nonrel P99.5-energy relative delta | PASS | 0.000702279 | 0.45 |
| rel P99.5-energy relative delta | PASS | 0.0346386 | 0.45 |
| nonrel P99.9-energy relative delta | PASS | 0.0102636 | 0.6 |
| rel P99.9-energy relative delta | PASS | 0.0245432 | 0.6 |
| nonrel fraction above 5Te relative delta | PASS | 0.00625903 | 0.35 |
| rel fraction above 5Te relative delta | PASS | 0.0115551 | 0.35 |
| nonrel fraction above 10Te relative delta | PASS | 0.00847458 | 0.75 |
| rel fraction above 10Te relative delta | PASS | 0.101695 | 0.75 |
| nonrel z-mean absolute delta | PASS | 0.0155901 | 0.25 |
| rel z-mean absolute delta | PASS | 0.0209603 | 0.25 |
| nonrel pitch-mean absolute delta | PASS | 0.000849295 | 0.05 |
| rel pitch-mean absolute delta | PASS | 8.8241e-05 | 0.05 |
| nonrel RF absorbed power ratio | PASS | 8.56214 | 25 |
| rel RF absorbed power ratio | PASS | 5.76517 | 25 |
| nonrel RF electric field positive | PASS | 0 | 0 |
| rel RF electric field positive | PASS | 0 | 0 |

Overall validation status: `PASS`.
