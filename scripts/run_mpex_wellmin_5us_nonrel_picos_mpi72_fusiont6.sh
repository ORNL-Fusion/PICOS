#!/usr/bin/env bash
set -euo pipefail

PICOS_ROOT=${PICOS_ROOT:-/home/78k/picos_kinetic_electron_eval/PICOS}
LINEAR_ROOT=${LINEAR_ROOT:-/home/78k/picos_kinetic_electron_eval/LinearFokkerPlanck_Axisymmetric}
MICROMAMBA=${MICROMAMBA:-/home/78k/picos_kinetic_electron_eval/micromamba/bin/micromamba}
PICOS_ENV=${PICOS_ENV:-/home/78k/picos_kinetic_electron_eval/envs/picos-ech}

RUN_LABEL=p16384_5us_nocoll_nonrel_mpi72_fusiont6
OUT_DIR=${OUT_DIR:-${PICOS_ROOT}/validation/mpex_scenario14_ex8_wellmin_65p588ghz_target8_${RUN_LABEL}}
PICOS_TAG=mpex_scenario14_ex8_well_min_65p588ghz_target8_${RUN_LABEL}_nonrel
CASE_NAME=xp_MPEX_Scenario14_Ex8_WellMin_65p588GHz_Target8_${RUN_LABEL}
DESCRIPTOR=ex8_well_min_65p588ghz_target8_${RUN_LABEL}

mkdir -p "${OUT_DIR}"
cd "${PICOS_ROOT}"

{
  date
  echo "Setting up PICOS++ nonrel MPEX scenario-14 well-minimum run"
  echo "Output directory: ${OUT_DIR}"
  echo "PICOS tag: ${PICOS_TAG}"

  env MPLBACKEND=Agg HDF5_USE_FILE_LOCKING=FALSE \
    "${MICROMAMBA}" run -p "${PICOS_ENV}" \
    python scripts/compare_mpex_scenario14_fig17.py \
      --setup \
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
      --particles 16384 \
      --steps 4000000 \
      --physical-time 5.0000000000000004e-06 \
      --collisions 0 \
      --omp-threads 1 \
      --fortran-omp-threads 48 \
      --picos-omp-threads 1 \
      --mpi-ranks 72 \
      --picos-run-mode nonrel \
      --case-name "${CASE_NAME}" \
      --descriptor "${DESCRIPTOR}" \
      --picos-tag-nonrel "${PICOS_TAG}" \
      --picos-tag-rel "${PICOS_TAG}_rel_unused" \
      --smooth-z-bins 8.0000000000000004e-01 \
      --smooth-energy-bins 8.0000000000000004e-01 \
      --mask-below-log-vmin \
      --out-dir "${OUT_DIR}"

  echo "Starting PICOS++ with 72 MPI ranks over hardware threads"
  cd "${PICOS_ROOT}/picosFILES"
  env HDF5_USE_FILE_LOCKING=FALSE OMP_NUM_THREADS=1 \
    "${MICROMAMBA}" run -p "${PICOS_ENV}" \
    mpirun --use-hwthread-cpus -np 72 \
      "${PICOS_ROOT}/build-fusiont6/picosFILES/src/xpicos" \
      1-D outputFiles "${PICOS_TAG}"

  cd "${PICOS_ROOT}"
  echo "PICOS++ finished; generating Fig. 17-style plot and summary"
  env MPLBACKEND=Agg HDF5_USE_FILE_LOCKING=FALSE \
    "${MICROMAMBA}" run -p "${PICOS_ENV}" \
    python scripts/compare_mpex_scenario14_fig17.py \
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
      --particles 16384 \
      --steps 4000000 \
      --physical-time 5.0000000000000004e-06 \
      --collisions 0 \
      --omp-threads 1 \
      --fortran-omp-threads 48 \
      --picos-omp-threads 1 \
      --mpi-ranks 72 \
      --picos-run-mode nonrel \
      --case-name "${CASE_NAME}" \
      --descriptor "${DESCRIPTOR}" \
      --picos-tag-nonrel "${PICOS_TAG}" \
      --picos-tag-rel "${PICOS_TAG}_rel_unused" \
      --smooth-z-bins 8.0000000000000004e-01 \
      --smooth-energy-bins 8.0000000000000004e-01 \
      --mask-below-log-vmin \
      --out-dir "${OUT_DIR}"
  date
} > "${OUT_DIR}/run.log" 2>&1
