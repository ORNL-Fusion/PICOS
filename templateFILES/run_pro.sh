#!/bin/bash

# Repository directory:
REPO_DIR=

# Directory for  libraries:
HDF5_INSTALL=/lib
ARMADILLO_INSTALL=/lib # /lib or /lib64

# Simulation ID
ID=""

# Dimensionality of simulation
DIMENSIONALITY="1-D"

# Available number of cores in system
NUM_CORES=32

# Number of MPI processes
NUM_MPI_PROCESSES=32

# Number of OMP threads per MPI
NUM_OMP_PER_MPI=$((NUM_CORES/NUM_MPI_PROCESSES))

# Location of outputs folder
LOC_OUTPUT_FOLDER=${REPO_DIR}"/picosFILES/outputFiles"

rm -r ${LOC_OUTPUT_FOLDER}"/"${ID}

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HDF5_INSTALL}:${ARMADILLO_INSTALL}
export HDF5_USE_FILE_LOCKING=FALSE

echo "LD_LIBRARY_PATH: "${LD_LIBRARY_PATH}
echo "Number of MPI processes: "${NUM_MPI_PROCESSES}
echo "Number of OMP threads per MPI: "${NUM_OMP_PER_MPI}


if [$ID == ""]; then
    echo "USING DEFAULT INPUT FILES"
    mpirun --use-hwthread-cpus -np $((NUM_MPI_PROCESSES)) -x LD_LIBRARY_PATH -x OMP_NUM_THREADS=$((NUM_OMP_PER_MPI)) -x HDF5_USE_FILE_LOCKING=FALSE bin/PICOS++ ${DIMENSIONALITY} ${LOC_OUTPUT_FOLDER}
else
    echo "USING MODIFIED INPUT FILES"
    mpirun --use-hwthread-cpus -np $((NUM_MPI_PROCESSES)) -x LD_LIBRARY_PATH -x OMP_NUM_THREADS=$((NUM_OMP_PER_MPI)) -x HDF5_USE_FILE_LOCKING=FALSE bin/PICOS++ ${DIMENSIONALITY} ${LOC_OUTPUT_FOLDER} 
${ID}
fi

cd   inputFiles/
cp input_file.input ions_properties.ion  *.txt ../outputFiles/
cd ..
cd outputFiles/
git log --oneline -1 > commitHash.txt
