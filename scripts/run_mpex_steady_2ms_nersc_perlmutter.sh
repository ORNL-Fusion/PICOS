#!/bin/bash
#SBATCH -A m77
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J mpex_steady2ms
#SBATCH --mail-user=kumara@ornl.gov
#SBATCH --mail-type=END,FAIL
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH -o slurm-%x-%j.out
#SBATCH -e slurm-%x-%j.err

set -euo pipefail

CASE_DIR=$(cd "$(dirname "$0")" && pwd)
PICOS_LOAD_MODULES=${PICOS_LOAD_MODULES:-1}
PICOS_GCC_MODULE=${PICOS_GCC_MODULE:-gcc-native/14}
if [ "${PICOS_LOAD_MODULES}" = "1" ] && command -v module >/dev/null 2>&1; then
  module load cpu
  module load PrgEnv-gnu
  module load "${PICOS_GCC_MODULE}"
  module load cray-hdf5
fi

if command -v g++ >/dev/null 2>&1; then
  LIBSTDCXX_DIR=$(dirname "$(g++ -print-file-name=libstdc++.so.6)")
  if [ -f "${LIBSTDCXX_DIR}/libstdc++.so.6" ]; then
    export LD_LIBRARY_PATH="${LIBSTDCXX_DIR}:${LD_LIBRARY_PATH:-}"
  fi
fi

PICOS_ROOT=${PICOS_ROOT:-${HOME}/myRepos/PICOS_ECH}
PICOS_BUILD_DIR=${PICOS_BUILD_DIR:-${PICOS_ROOT}/build}
PICOS_BIN=${PICOS_BIN:-${PICOS_BUILD_DIR}/picosFILES/src/xpicos}
RUN_ROOT=${RUN_ROOT:-${SCRATCH}/picosRuns/MPEX_ECH_runs/PICOS_NERSC_MPEX_steady_restart_2ms}
CASE_TAG=${CASE_TAG:-mpex_scenario14_ex8_steady_2ms_p262144_coll_mpexprof_nersc_nonrel}
MPI_RANKS=${MPI_RANKS:-128}
CPUS_PER_TASK=${CPUS_PER_TASK:-1}

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
  exit 2
fi
if [ $((MPI_RANKS % 2)) -ne 0 ]; then
  echo "MPI_RANKS must be even for PICOS++." >&2
  exit 2
fi

RUN_PICOS_FILES=${RUN_ROOT}/picosFILES
mkdir -p "${RUN_PICOS_FILES}/inputFiles" "${RUN_PICOS_FILES}/outputFiles" "${RUN_ROOT}/logs"
if [ -d "${PICOS_ROOT}/templateFILES" ]; then
  rsync -a "${PICOS_ROOT}/templateFILES/" "${RUN_PICOS_FILES}/inputFiles/"
fi
if [ -d "${PICOS_ROOT}/picosFILES/inputFiles" ]; then
  rsync -a "${PICOS_ROOT}/picosFILES/inputFiles/" "${RUN_PICOS_FILES}/inputFiles/"
fi
if [ -d "${CASE_DIR}/inputFiles" ]; then
  rsync -a "${CASE_DIR}/inputFiles/" "${RUN_PICOS_FILES}/inputFiles/"
fi
git -C "${PICOS_ROOT}" log --oneline -1 > "${RUN_ROOT}/commitHash.txt"

input_file="${RUN_PICOS_FILES}/inputFiles/input_file_${CASE_TAG}.input"
ion_file="${RUN_PICOS_FILES}/inputFiles/ions_properties_${CASE_TAG}.ion"
if [ ! -f "${input_file}" ] || [ ! -f "${ion_file}" ]; then
  echo "Missing input deck for ${CASE_TAG}." >&2
  echo "Regenerate it with scripts/compare_mpex_scenario14_fig17.py --setup --rf-heating 0 --physical-time 2.0e-3." >&2
  exit 2
fi

if [ -d "${RUN_PICOS_FILES}/outputFiles/${CASE_TAG}" ]; then
  if [ "${OVERWRITE_OUTPUTS:-0}" = "1" ]; then
    rm -rf "${RUN_PICOS_FILES}/outputFiles/${CASE_TAG}"
  else
    echo "Output directory already exists for ${CASE_TAG}; set OVERWRITE_OUTPUTS=1 to replace it." >&2
    exit 2
  fi
fi

echo "[$(date)] Starting ${CASE_TAG}"
echo "PICOS_ROOT=${PICOS_ROOT}"
echo "PICOS_BIN=${PICOS_BIN}"
echo "RUN_ROOT=${RUN_ROOT}"
echo "RUN_PICOS_FILES=${RUN_PICOS_FILES}"
echo "MPI_RANKS=${MPI_RANKS}"
echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "LIBSTDCXX_DIR=${LIBSTDCXX_DIR:-unset}"

(
  cd "${RUN_PICOS_FILES}"
  srun -n "${MPI_RANKS}" -c "${CPUS_PER_TASK}" --cpu-bind=cores \
    "${PICOS_BIN}" 1-D outputFiles "${CASE_TAG}"
) > "${RUN_ROOT}/logs/${CASE_TAG}.log" 2>&1

restart_path="${RUN_PICOS_FILES}/outputFiles/${CASE_TAG}/HDF5"
echo "[$(date)] Finished ${CASE_TAG}"
echo "Restart path for ECH continuation:"
echo "${restart_path}"
