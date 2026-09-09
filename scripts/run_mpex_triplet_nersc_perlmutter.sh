#!/bin/bash
#SBATCH -A m77
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J mpex_ech_triplet
#SBATCH --mail-user=kumara@ornl.gov
#SBATCH --mail-type=END,FAIL
#SBATCH -t 11:30:00
#SBATCH -o slurm-%x-%j.out
#SBATCH -e slurm-%x-%j.err

set -euo pipefail

if [ "${PE_ENV:-}" == "CRAY" ]; then
  export FILENV=my_filenenv
  assign -U on g:sf
fi

# Source tree is expected in myRepos; run products go to scratch.
PICOS_ROOT=${PICOS_ROOT:-${HOME}/myRepos/PICOS}
PICOS_BUILD_DIR=${PICOS_BUILD_DIR:-${PICOS_ROOT}/build}
PICOS_BIN=${PICOS_BIN:-${PICOS_BUILD_DIR}/picosFILES/src/xpicos}
RUN_ROOT=${RUN_ROOT:-${SCRATCH}/PICOS_MPEX/scenario14_triplet}
RUN_LABEL=${RUN_LABEL:-p262144_100us_nocoll_nersc}

# Perlmutter CPU defaults. Keep MPI ranks even; PICOS++ requires that.
MPI_RANKS=${MPI_RANKS:-128}
CPUS_PER_TASK=${CPUS_PER_TASK:-2}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OMP_PLACES=${OMP_PLACES:-threads}
export OMP_PROC_BIND=${OMP_PROC_BIND:-spread}
export HDF5_USE_FILE_LOCKING=FALSE

HDF5_INSTALL=${HDF5_INSTALL:-${PICOS_ROOT}/HDF5/lib}
ARMADILLO_INSTALL=${ARMADILLO_INSTALL:-${PICOS_ROOT}/arma_libs/lib64}
if [ -d "${HDF5_INSTALL}" ]; then
  export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:${HDF5_INSTALL}"
fi
if [ -d "${ARMADILLO_INSTALL}" ]; then
  export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:${ARMADILLO_INSTALL}"
fi

if [ ! -x "${PICOS_BIN}" ]; then
  echo "Missing PICOS++ executable: ${PICOS_BIN}" >&2
  echo "Set PICOS_BUILD_DIR or PICOS_BIN to the build under ${PICOS_ROOT}." >&2
  exit 2
fi

if [ $((MPI_RANKS % 2)) -ne 0 ]; then
  echo "MPI_RANKS must be even for PICOS++." >&2
  exit 2
fi

RUN_PICOS_FILES=${RUN_ROOT}/picosFILES
mkdir -p "${RUN_PICOS_FILES}/inputFiles" "${RUN_PICOS_FILES}/outputFiles" "${RUN_ROOT}/logs"
rsync -a "${PICOS_ROOT}/picosFILES/inputFiles/" "${RUN_PICOS_FILES}/inputFiles/"
git -C "${PICOS_ROOT}" log --oneline -1 > "${RUN_ROOT}/commitHash.txt"

CASES=(
  "mpex_scenario14_ex8_left_resonance_70ghz_target8_${RUN_LABEL}_nonrel"
  "mpex_scenario14_ex8_right_resonance_70ghz_target8_${RUN_LABEL}_nonrel"
  "mpex_scenario14_ex8_well_min_65p588ghz_target8_${RUN_LABEL}_nonrel"
)

echo "PICOS_ROOT=${PICOS_ROOT}"
echo "PICOS_BIN=${PICOS_BIN}"
echo "RUN_ROOT=${RUN_ROOT}"
echo "RUN_LABEL=${RUN_LABEL}"
echo "MPI_RANKS=${MPI_RANKS}"
echo "CPUS_PER_TASK=${CPUS_PER_TASK}"
echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "Cases: ${CASES[*]}"

for tag in "${CASES[@]}"; do
  input_file="${RUN_PICOS_FILES}/inputFiles/input_file_${tag}.input"
  ion_file="${RUN_PICOS_FILES}/inputFiles/ions_properties_${tag}.ion"
  if [ ! -f "${input_file}" ] || [ ! -f "${ion_file}" ]; then
    echo "Missing input deck for ${tag}." >&2
    echo "Generate decks before submitting, for example:" >&2
    echo "  python scripts/run_mpex_triplet_picos_fortran.py --setup-only --skip-smooth --picos-run-mode nonrel --run-label ${RUN_LABEL} --particles 262144 --physical-time 1.0e-4 --output-count 50 --collisions 0" >&2
    exit 2
  fi

  if [ -d "${RUN_PICOS_FILES}/outputFiles/${tag}" ]; then
    if [ "${OVERWRITE_OUTPUTS:-0}" == "1" ]; then
      rm -rf "${RUN_PICOS_FILES}/outputFiles/${tag}"
    else
      echo "Output directory already exists for ${tag}; set OVERWRITE_OUTPUTS=1 to replace it." >&2
      exit 2
    fi
  fi

  echo "[$(date)] Starting ${tag}"
  (
    cd "${RUN_PICOS_FILES}"
    srun -n "${MPI_RANKS}" -c "${CPUS_PER_TASK}" --cpu-bind=cores \
      "${PICOS_BIN}" 1-D outputFiles "${tag}"
  ) > "${RUN_ROOT}/logs/${tag}.log" 2>&1
  echo "[$(date)] Finished ${tag}"
done

echo "All MPEX scenario-14 triplet PICOS++ runs completed."
echo "Outputs: ${RUN_PICOS_FILES}/outputFiles"
