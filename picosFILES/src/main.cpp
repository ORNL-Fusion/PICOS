// PICOS++
// Include intrinsic header files:
// =============================================================================
#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>
#include <ctime>
#include <typeinfo>
#include <utility>

// Include user-defined header files:
// =============================================================================
#include "types.h"
/*
#include "collisionOperator.h"
#include "particleBoundaryConditions.h"
#include "rfOperator.h"
#include "initialize.h"
#include "PIC.h"
#include "fields.h"
#include "units.h"
#include "outputHDF5.h"
*/

// Include headers for parallelization:
// =============================================================================
#include <omp.h>
//#include "mpi_main.h"

using namespace std;
using namespace arma;

// #############################################################################
// MAIN function:
// #############################################################################

int main(int argc, char* argv[])
{
    /*
    MPI_Init(&argc, &argv);
    */

    // Create simulation objects:
    // =========================================================================
    // MPI object to hold topology information:
    /*
    MPI_MAIN mpi_main;
    */

    // Input parameters for simulation:
    /*
    simulationParameters params;
    */

    // Ion species vector:
    vector<ionSpecies_TYP> IONS;

    cout << "Hello from PICOS++!!" << endl;

	return(0);
}
