#!/bin/bash

REPO_DIR=${PWD}

# * * * * * Please change these variables as required * * * * *
export CC=gcc
export CPP=cpp
export CXX=g++
export FC=gfortran

MPICXX=mpic++

HDF5_INSTALLATION_FOLDER=$REPO_DIR
ARMADILLO_INSTALLATION_FOLDER=$REPO_DIR
# * * * * * Please change these variables as required * * * * *

HDF5_INSTALLATION_FOLDER=$HDF5_INSTALLATION_FOLDER"/HDF5"
ARMADILLO_INSTALLATION_FOLDER=$ARMADILLO_INSTALLATION_FOLDER"/arma_libs"

HDF5_VERSION='hdf5-1.10.4'
ARMADILLO_VERSION='armadillo-9.850.1'

# Delete any existing previous instalation of the libraries
rm -r $HDF5_INSTALLATION_FOLDER
rm -r $ARMADILLO_INSTALLATION_FOLDER

# Create new folders for local installation of HDF5 and armadillo
mkdir $HDF5_INSTALLATION_FOLDER
mkdir $ARMADILLO_INSTALLATION_FOLDER

# Delete existing HDF5 and armadillo build folders
rm -r $HDF5_VERSION
rm -r $ARMADILLO_VERSION

# Local installation of Armadillo library

# We extract the source code of armadillo library
tar -xvf $ARMADILLO_VERSION$".tar.xz"

# We enter the source code directory of armadillo
cd $ARMADILLO_VERSION

# Armadillo installation
cmake -DCMAKE_INSTALL_PREFIX=$ARMADILLO_INSTALLATION_FOLDER -DDETECT_HDF5=false .

make

if [ $? -eq 0 ] ; then
make install

cd ../

rm -r $ARMADILLO_VERSION
else
echo 'ERROR: Uh-oh! Something went wrong in ARMADILLO installation!'
return
fi

# Local installation of HDF5 library

# We extract the HDF5 source code
tar -xvf $HDF5_VERSION".tar.gz"

# Enter the HDF5 source code folder
cd $HDF5_VERSION

# Set the prefix for installation folder
PREFIX=$HDF5_INSTALLATION_FOLDER

./configure --prefix=$PREFIX --enable-cxx --enable-build-mode=production

make

if [ $? -eq 0 ] ; then
make install

cd ../

rm -r $HDF5_VERSION
else
echo 'ERROR: Uh-oh! Something went wrong in HDF5 installation!'
return
fi

# Setting up Makefile and Compile environment variables

rm picosFILES/compile.sh
rm picosFILES/Makefile
rm picosFILES/run_pro.sh
rm picosFILES/inputFiles/input_file.input
rm picosFILES/inputFiles/ions_properties.ion

mkdir picosFILES/inputFiles
mkdir picosFILES/outputFiles
mkdir picosFILES/bin
mkdir picosFILES/obj

cp templateFILES/compile.sh picosFILES/compile.sh
cp templateFILES/Makefile picosFILES/Makefile
cp templateFILES/run_pro.sh picosFILES/run_pro.sh
cp templateFILES/debug_run_pro.sh picosFILES/debug_run_pro.sh
cp templateFILES/profile_run_pro.sh picosFILES/profile_run_pro.sh
cp templateFILES/input_file.input picosFILES/inputFiles/input_file.input
cp templateFILES/ions_properties.ion picosFILES/inputFiles/ions_properties.ion


MPICXX_=$(echo ${MPICXX} | sed -e 's/\//_SEP_/g')
sed -i 's/MPICXX=/'"MPICXX=${MPICXX_}/g" picosFILES/Makefile

HDF5_INSTALLATION_FOLDER_=$(echo ${HDF5_INSTALLATION_FOLDER} | sed -e 's/\//_SEP_/g')
sed -i 's/HDF5_INSTALL=/'"HDF5_INSTALL=${HDF5_INSTALLATION_FOLDER_}/g" picosFILES/Makefile
sed -i 's/HDF5_INSTALL=/'"HDF5_INSTALL=${HDF5_INSTALLATION_FOLDER_}/g" picosFILES/run_pro.sh
sed -i 's/HDF5_INSTALL=/'"HDF5_INSTALL=${HDF5_INSTALLATION_FOLDER_}/g" picosFILES/debug_run_pro.sh
sed -i 's/HDF5_INSTALL=/'"HDF5_INSTALL=${HDF5_INSTALLATION_FOLDER_}/g" picosFILES/profile_run_pro.sh

ARMADILLO_INSTALLATION_FOLDER_=$(echo ${ARMADILLO_INSTALLATION_FOLDER} | sed -e 's/\//_SEP_/g')
sed -i 's/ARMADILLO_INSTALL=/'"ARMADILLO_INSTALL=${ARMADILLO_INSTALLATION_FOLDER_}/g" picosFILES/Makefile
sed -i 's/ARMADILLO_INSTALL=/'"ARMADILLO_INSTALL=${ARMADILLO_INSTALLATION_FOLDER_}/g" picosFILES/run_pro.sh
sed -i 's/ARMADILLO_INSTALL=/'"ARMADILLO_INSTALL=${ARMADILLO_INSTALLATION_FOLDER_}/g" picosFILES/debug_run_pro.sh
sed -i 's/ARMADILLO_INSTALL=/'"ARMADILLO_INSTALL=${ARMADILLO_INSTALLATION_FOLDER_}/g" picosFILES/profile_run_pro.sh

REPO_DIR_=$(echo ${REPO_DIR} | sed -e 's/\//_SEP_/g')
sed -i 's/REPO_DIR=/'"REPO_DIR=${REPO_DIR_}/g" picosFILES/run_pro.sh
sed -i 's/REPO_DIR=/'"REPO_DIR=${REPO_DIR_}/g" picosFILES/debug_run_pro.sh
sed -i 's/REPO_DIR=/'"REPO_DIR=${REPO_DIR_}/g" picosFILES/profile_run_pro.sh

sed -i 's/_SEP_/\//g' picosFILES/Makefile
sed -i 's/_SEP_/\//g' picosFILES/run_pro.sh
sed -i 's/_SEP_/\//g' picosFILES/debug_run_pro.sh
sed -i 's/_SEP_/\//g' picosFILES/profile_run_pro.sh

sed -i 's/ARMA_LIBS/'"${ARMADILLO_INSTALLATION_FOLDER_}/g" picosFILES/compile.sh
sed -i 's/_SEP_/\//g' picosFILES/compile.sh

# Compile for the first times
cd picosFILES

./compile.sh

if [ $? -eq 0 ] ; then
echo ''
echo '* * * * * * * * * * * * * * * * * * * * * *'
echo '*          Installation succeeded         *'
echo '* * * * * * * * * * * * * * * * * * * * * *'
echo ''
fi
