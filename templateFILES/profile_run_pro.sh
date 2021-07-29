#!/bin/bash

REPO_DIR=
HDF5_INSTALL=/lib
ARMADILLO_INSTALL=/lib
FORGE_PATH=

# Simulation ID
ID=""

# Dimensionality of simulation
DIMENSIONALITY="1-D"

# Available number of cores in system
NUM_CORES=56

# Number of MPI processes
NUM_MPI_PROCESSES=14

# Number of OMP threads per MPI
NUM_OMP_PER_MPI=$((NUM_CORES/NUM_MPI_PROCESSES))

# Location of outputs folder
LOC_OUTPUT_FOLDER=${REPO_DIR}"/picosFILES/outputFiles"

rm -r ${LOC_OUTPUT_FOLDER}"/"${ID}

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HDF5_INSTALL}:${ARMADILLO_INSTALL}
export PATH=${FORGE_PATH}:$PATH

echo "LD_LIBRARY_PATH: "${LD_LIBRARY_PATH}
echo "Number of MPI processes: "${NUM_MPI_PROCESSES}
echo "Number of OMP threads per MPI: "${NUM_OMP_PER_MPI}


if [$ID == ""]; then
    echo "USING DEFAULT INPUT FILES"
    map --profile mpirun --use-hwthread-cpus -np $((NUM_MPI_PROCESSES)) -x LD_LIBRARY_PATH -x OMP_NUM_THREADS=$((NUM_OMP_PER_MPI)) bin/PICOS++ ${DIMENSIONALITY} ${LOC_OUTPUT_FOLDER}
else
    echo "USING MODIFIED INPUT FILES"
    map --profile mpirun --use-hwthread-cpus -np $((NUM_MPI_PROCESSES)) -x LD_LIBRARY_PATH -x OMP_NUM_THREADS=$((NUM_OMP_PER_MPI)) bin/PICOS++ ${DIMENSIONALITY} ${LOC_OUTPUT_FOLDER} ${ID}
fi
