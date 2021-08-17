#Makefile
UNAME:=$(shell uname)


SRC=src
OBJ=obj
LIBS=libs
BIN=bin

MPICXX=
CCXXFLAGS=-std=c++11
CCDEFINEFLAGS=-DHDF5_FLOAT	-DHDF5_SINGLE_FILE	-DCHECKS_ON
CCFLAGS=$(CCXXFLAGS)	-fopenmp	-O3
HDF5_INSTALL=
ARMADILLO_INSTALL=

INCLUDE_HDF5=-I$(HDF5_INSTALL)/include
INCLUDES=-I$(ARMADILLO_INSTALL)/include	-I$(ARMADILLO_INSTALL)/lib	-I$(ARMADILLO_INSTALL)/lib64
EXTLIB=-L$(HDF5_INSTALL)/lib
LIBSHDF=$(EXTLIB) $(HDF5_INSTALL)/lib/libhdf5.a
LIB_FLAGS=-L$(ARMADILLO_INSTALL)/lib	-L$(ARMADILLO_INSTALL)/lib64
LIB=-lz	-ldl	-lm	-lhdf5_cpp	-larmadillo #-lgfortran
OMP_LIB_FLAGS=-lgomp

BINARIES = $(addprefix $(OBJ)/,\
main.o	            \
types.o   			    \
initialize.o		    \
units.o						  \
initDistribution.o  \
outputHDF5.o 			  \
PIC.o 						  \
particleBC.o			  \
collisionOperator.o \
mpi_main.o)

all: $(BIN)/PICOS++

$(OBJ)/%.o : $(SRC)/%.cpp
	$(MPICXX)	-c	$(CCFLAGS)	$(CCDEFINEFLAGS)	$<	$(INCLUDES)	$(INCLUDE_HDF5)	$(LIB_FLAGS) -o $@

$(BIN)/PICOS++: $(BINARIES)
	$(MPICXX)	$(BINARIES)	-o	$@	$(LIB_FLAGS)	$(LIBSHDF)	$(LIB)	$(OMP_LIB_FLAGS)

info:
	@echo	"COMPILING IN: "$(UNAME)

clean:
	@echo	"CLEANING PREVIOUS COMPILATION"
	rm -f $(BIN)/PICOS++	$(OBJ)/*.o	$(SRC)/*~
