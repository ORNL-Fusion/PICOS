#!/usr/bin/env bash
set -euo pipefail

PICOS_ROOT=${PICOS_ROOT:-/home/78k/picos_kinetic_electron_eval/PICOS}
LINEAR_ROOT=${LINEAR_ROOT:-/home/78k/picos_kinetic_electron_eval/LinearFokkerPlanck_Axisymmetric}
MICROMAMBA=${MICROMAMBA:-/home/78k/picos_kinetic_electron_eval/micromamba/bin/micromamba}
PICOS_ENV=${PICOS_ENV:-/home/78k/picos_kinetic_electron_eval/envs/picos-ech}

OUT_DIR=${OUT_DIR:-${PICOS_ROOT}/validation/mpex_scenario14_ex8_wellmin_65p588ghz_target8_p2048_1us_nocoll_nonrel_mpi48_fusiont6}

# Reuse the completed Fortran output from the first 1 us attempt.  Only the
# PICOS++ nonrel stage failed there because OpenMPI rejected 64 requested ranks.
FORTRAN_CASE=xp_MPEX_Scenario14_Ex8_WellMin_65p588GHz_Target8_p2048_1us_nocoll_nonrel_mpi64_fusiont6
FORTRAN_DESCRIPTOR=ex8_well_min_65p588ghz_target8_p2048_1us_nocoll_nonrel_mpi64_fusiont6
PICOS_TAG=mpex_scenario14_ex8_well_min_65p588ghz_target8_p2048_1us_nocoll_nonrel_mpi48_fusiont6_nonrel

mkdir -p "${OUT_DIR}"
cd "${PICOS_ROOT}"

exec env MPLBACKEND=Agg HDF5_USE_FILE_LOCKING=FALSE \
  "${MICROMAMBA}" run -p "${PICOS_ENV}" \
  python scripts/compare_mpex_scenario14_fig17.py \
    --setup \
    --run-picos \
    --plot \
    --linear-root "${LINEAR_ROOT}" \
    --picos-build-dir "${PICOS_ROOT}/build-fusiont6" \
    --scenario14-b-file "${PICOS_ROOT}/templateFILES/MPEX_B_norm_PICOS_scenario_14.txt" \
    --b-scale 2.0041865242724399e-01 \
    --z-min=-2.0000000000000000e+00 \
    --z-max=8.0000000000000000e+00 \
    --plot-z-min=-2.0000000000000000e+00 \
    --plot-z-max=8.0000000000000000e+00 \
    --energy-max 3.5000000000000000e+03 \
    --rf-frequency 6.5587790285535164e+10 \
    --no-derive-rf-frequency-from-resonance-b \
    --resonance-b-t 1.1715247658329999e+00 \
    --resonance-z 2.8743718590000000e+00 \
    --rf-window-end 8.0000000000000000e+00 \
    --target-z 8.0000000000000000e+00 \
    --particles 2048 \
    --steps 800000 \
    --physical-time 1.0000000000000000e-06 \
    --collisions 0 \
    --omp-threads 1 \
    --fortran-omp-threads 48 \
    --picos-omp-threads 1 \
    --mpi-ranks 48 \
    --picos-run-mode nonrel \
    --case-name "${FORTRAN_CASE}" \
    --descriptor "${FORTRAN_DESCRIPTOR}" \
    --picos-tag-nonrel "${PICOS_TAG}" \
    --picos-tag-rel "${PICOS_TAG}_rel_unused" \
    --smooth-z-bins 8.0000000000000004e-01 \
    --smooth-energy-bins 8.0000000000000004e-01 \
    --mask-below-log-vmin \
    --out-dir "${OUT_DIR}" \
  > "${OUT_DIR}/run.log" 2>&1
