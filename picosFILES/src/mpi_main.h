#ifndef H_MPI_MAIN
#define H_MPI_MAIN

#include <iostream>
#include <vector>
#include <armadillo>
#include <omp.h>
#include <cmath>

#include "types.h"

#include "mpi.h"

using namespace std;
using namespace arma;

class MPI_MAIN_TYP
{

private:


public:

	//! Constructor of MPI_MAIN
	MPI_MAIN_TYP(){};

	void createMPITopology(params_TYP * params);

	void finalizeCommunications(params_TYP * params);

};

#endif
