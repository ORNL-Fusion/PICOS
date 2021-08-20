// PICOS++
// Include intrinsic header files:
// =============================================================================
#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>
#include <ctime>
#include <utility>

// Include user-defined header files:
// =============================================================================
#include "types.h"
#include "initialize.h"
#include "units.h"
#include "outputHDF5.h"
#include "PIC.h"
#include "particleBC.h"
#include "collisionOperator.h"
#include "fieldSolve.h"
#include "rfOperator.h"

// Include headers for parallelization:
// =============================================================================
#include <omp.h>
#include "mpi_main.h"

using namespace std;
using namespace arma;

// #############################################################################
// MAIN function:
// #############################################################################

int main(int argc, char* argv[])
{
    // Initialize MPI process:
    // =========================================================================
    MPI_Init(&argc, &argv);

    // Create simulation objects:
    // =========================================================================
    // MPI object to hold topology information:
    MPI_MAIN_TYP mpi_main;

    // Input parameters for simulation:
    params_TYP params;

    // Ion species vector:
    vector<ionSpecies_TYP> IONS;

    // Characteristic scales:
    CS_TYP CS;

    // Electromagnetic fields:
    fields_TYP fields;

    // Collision operator object:
    coll_operator_TYP coll_operator;

    // Particle boundary condition operator:
    particleBC_TYP particleBC;

    // Initialize "params" based on input file:
    init_TYP init(&params, argc, argv);

    // UNITS object of type "FT" and "FT":
    units_TYP units;

    // Initialize simulation objects:
    // =========================================================================
    // Create MPI topology:
    mpi_main.createMPITopology(&params);

    // Read "ions_properties.ion" and populate "IONS" vector
    init.loadIonParameters(&params, &IONS);

    // Define characteristic scales and broadcast them to all processes in COMM_WORLD:
    units.defineCharacteristicScalesAndBcast(&params, &IONS, &CS);

    // Create object and allocate memory according to "params":
    FS_TYP FS(&params);

    // Define fundamental scales and broadcast them to all processes in COMM_WORLD:
    units.calculateFundamentalScalesAndBcast(&params, &IONS, &FS);

    // Define mesh geometry and populate "params" (FS is not used):
    init.loadMeshGeometry(&params, &FS);

    // Read IC external profile and load them to "params.p_IC":
    init.loadPlasmaProfiles(&params, &IONS);

    // Check that mesh size is consistent with hybrid approximation:
    units.spatialScalesSanityCheck(&params, &FS);

    // Initialize electromagnetic field variable:
    init.initializeFields(&params, &fields);

    // Initialize IONS: scalar, bulk and particle arrays
    init.setupIonsInitialCondition(&params, &CS, &fields, &IONS);

    // HDF object constructor and create "main.h5"
    HDF_TYP HDF(&params, &FS, &IONS);

    // Define time step based on ion CFL condition:
    units.defineTimeStep(&params, &IONS);
    //cout << "checkpoint 3" << endl;

    // Normalize "params", "IONS", "fields" using "CS"
    units.normalizeVariables(&params, &IONS, &fields, &CS);

    // #########################################################################
    /**************** All the quantities below are dimensionless ****************/
    // #########################################################################

    // Definition of variables for advancing in time particles and fields:
    // =========================================================================
    double t1 = 0.0;
    double t2 = 0.0;
    params.currentTime = 0.0;
    int outputIterator = 0;
    int numberOfIterationsForEstimator = 1000;

    // Create EM solver:
    // =========================================================================
    fields_solver_TYP fields_solver(&params, &CS); // Initializing the EMF_SOLVER class object.

    // Create PIC solver:
    // =========================================================================
    PIC_TYP PIC(&params, &CS, &fields, &IONS);
    /*
    // Run 3 dummy cycles to load "n" and "nv" at previous time steps:
    // =========================================================================
    for(int tt=0; tt<3; tt++)
    {
        PIC.extrapolateIonsMoments(&params, &fields, &IONS);
        PIC.extrapolateMoments_AllSpecies(&params, &fields, &IONS);
    }
    */

    // Create RF operator object:
    // =========================================================================
    RF_Operator_TYP RF_operator(&params,&CS,&fields,&IONS);

    // Save 1st output:
    // =========================================================================
    HDF.saveOutputs(&params, &IONS, &fields, &CS, 0, 0);

    // Start timing simulations:
    // =========================================================================
    t1 = MPI_Wtime();

    // #########################################################################
    // Start time iterations:
    // #########################################################################
    for(int tt=0; tt<params.timeIterations; tt++)
    {

        if (params.mpi.IS_PARTICLES_ROOT)
        {
            if (fmod((double)(tt + 1), 1000) == 0)
            {
                cout << "time = " << tt*params.DT*CS.time*1E3 << " [ms] "<< endl;

            }
        }

        // Advance particles and re-inject:
        // =====================================================================
        if (params.SW.advancePos == 1)
        {
            // Advance particle position and velocity to level X^(N+1):
            PIC.advanceParticles(&params, &fields, &IONS);

            // Re-inject particles that leave computational domain:
            particleBC.applyParticleReinjection(&params,&CS,&fields,&IONS);

            // Assign cell:
            PIC.assignCell_AllSpecies(&params,&IONS);

            // Interpolate all fields:
            PIC.interpolateFields_AllSpecies(&params,&IONS,&fields);
        }

        // Calculate ion moments:
        // =====================================================================
        PIC.extrapolateMoments_AllSpecies(&params,&CS,&fields,&IONS);

        // Apply collision operator:
        // =====================================================================
        if (params.SW.Collisions == 1)
        {
            coll_operator.ApplyCollisions_AllSpecies(&params,&CS,&IONS);
        }

        // Apply RF operator:
        // =====================================================================
        if (params.SW.RFheating == 1)
        {
            RF_operator.ApplyRfHeating(&params,&CS,&fields,&IONS);
        }

        // Field solve:
        // =====================================================================
        // Magnetic field:
        if (params.SW.BfieldSolve == 1)
        {
            // Use Faraday's law to advance the magnetic field to level B^(N+1).
            //fields_solver.advanceBfield(&params, &fields, &IONS);
        }

        if (params.SW.EfieldSolve == 1)
        {
            // Use Ohm's law to advance the electric field:
            fields_solver.advanceEfield(&params,&fields,&CS,&IONS);
        }

        // Advance time:
        // =====================================================================
        params.currentTime += params.DT*CS.time;

        // Save data:
        // =====================================================================
        if(fmod((double)(tt + 1), params.outputCadenceIterations) == 0)
        {
            vector<ionSpecies_TYP> IONS_OUT = IONS;

            HDF.saveOutputs(&params, &IONS_OUT, &fields, &CS, outputIterator+1, params.currentTime);

            outputIterator++;
        }

        // Estimate simulation time:
        // =====================================================================
        if(tt == numberOfIterationsForEstimator)
        {
            t2 = MPI_Wtime();

            double estimatedSimulationTime = ( (double)params.timeIterations*(t2 - t1)/(double)numberOfIterationsForEstimator )/60.0;

            if(params.mpi.MPI_DOMAIN_NUMBER == 0)
            {
                cout << "ESTIMATED TIME OF COMPLETION: " << estimatedSimulationTime <<" MINUTES" << endl;
            }
        }

    }
    // #########################################################################
    // End time iterations.
    // #########################################################################

    // Finalizing MPI communications:
    // =========================================================================
    mpi_main.finalizeCommunications(&params);

    cout << "end" << endl;

	return(0);
}
