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
#include "mpi_main.h"

using namespace std;
using namespace arma;

// #############################################################################
// MAIN function:
// #############################################################################

int main(int argc, char* argv[])
{
    // Initialize MPI process:
    MPI_Init(&argc, &argv);

    // Create simulation objects:
    // =========================================================================
    // MPI object to hold topology information:
    MPI_MAIN_TYP mpi_main;

    // Input parameters for simulation:
    params_TYP params;

    // temp stuff:
    params.mpi.MPIS_FIELDS = 2;
    params.timeIterations = 100000;
    params.SW.advancePos == 1;
    params.SW.Collisions == 1;
    params.SW.BfieldSolve == 1;
    params.SW.EfieldSolve == 1;
    MPI_Comm_size(MPI_COMM_WORLD, &params.mpi.NUMBER_MPI_DOMAINS);
    MPI_Comm_rank(MPI_COMM_WORLD, &params.mpi.MPI_DOMAIN_NUMBER);

    // Ion species vector:
    vector<ionSpecies_TYP> IONS;

    // Characteristic scales:
    CS_TYP CS;

    // Electromagnetic fields:
    fields_TYP fields;

    // Collision operator object:
    //collisionOperator FPCOLL;

    // Particle boundary condition operator:
    //PARTICLE_BC particleBC;

    // Initialize "params" based on input file:
    //INITIALIZE<IT, FT> init(&params, argc, argv);

    // UNITS object of type "FT" and "FT":
    //UNITS<IT, FT> units;

    // Initialize simulation objects:
    // =========================================================================
    // Create MPI topology:
    mpi_main.createMPITopology(&params);

    // Read "ions_properties.ion" and populate "IONS" vector
    //init.loadIonParameters(&params, &IONS);
    // Define characteristic scales and broadcast them to all processes in COMM_WORLD:
    //units.defineCharacteristicScalesAndBcast(&params, &IONS, &CS);
    // Create object and allocate memory according to "params":
    //fundamentalScales FS(&params);
    // Define fundamental scales and broadcast them to all processes in COMM_WORLD:
    //units.calculateFundamentalScalesAndBcast(&params, &IONS, &FS);
    // Define mesh geometry and populate "params" (FS is not used):
    //init.loadMeshGeometry(&params, &FS);
    // Read external profile files and load them to "params":
    //init.loadPlasmaProfiles(&params, &IONS);
    // Check that mesh size is consistent with hybrid approximation:
    //units.spatialScalesSanityCheck(&params, &FS);
    // Initialize electromagnetic field variable:
    //init.initializeFields(&params, &EB);
    // Initialize IONS: scalar, bulk and particle arrays
    //init.setupIonsInitialCondition(&params, &CS, &EB, &IONS);
    // HDF object constructor:
    //HDF<IT, FT> hdfObj(&params, &FS, &IONS);
    // Define time step based on CFL condition: Whistler and ion velocity
    //units.defineTimeStep(&params, &IONS);
    // Normalize "params", "IONS", "EB" using "CS"
    //units.normalizeVariables(&params, &IONS, &EB, &CS);

    // #########################################################################
    /**************** All the quantities below are dimensionless ****************/
    // #########################################################################

    // Definition of variables for advancing in time particles and fields:
    // =========================================================================
    double t1 = 0.0;
    double t2 = 0.0;
    double currentTime = 0.0;
    int outputIterator = 0;
    int numberOfIterationsForEstimator = 1000;

    // Create objects:
    // =========================================================================
    //EMF_SOLVER fields_solver(&params, &CS); // Initializing the EMF_SOLVER class object.
    //PIC ionsDynamics; // Initializing the PIC class object.

    // Run 3 dummy cycles to load "n" and "nv" at previous time steps:
    // =========================================================================
    for(int tt=0; tt<3; tt++)
    {
        //ionsDynamics.advanceIonsPosition(&params, &EB, &IONS, 0);

        //ionsDynamics.advanceIonsVelocity(&params, &CS, &EB, &IONS, 0);

        //ionsDynamics.extrapolateIonsMoments(&params, &EB, &IONS);
    }

    // Save 1st output:
    // =========================================================================
    //hdfObj.saveOutputs(&params, &IONS, &EB, &CS, 0, 0);

    // Start timing simulations:
    // =========================================================================
    t1 = MPI_Wtime();


    // #########################################################################
    // Start time iterations:
    // #########################################################################
    for(int tt=0; tt<params.timeIterations; tt++)
    {
        // Advance velocity:
        // =====================================================================
        if(tt == 0)
        {
             // Initial condition time level V^(1/2):
            //ionsDynamics.advanceIonsVelocity(&params, &CS, &EB, &IONS, 0.5*params.DT);
        }
        else
        {
             // Advance ions' velocity V^(N+1/2):
            //ionsDynamics.advanceIonsVelocity(&params, &CS, &EB, &IONS, params.DT);
        }

        // Advance position:
        // =====================================================================
        if (params.SW.advancePos == 1)
        {
            // Advance ions' position in time to level X^(N+1).
            //ionsDynamics.advanceIonsPosition(&params,&EB, &IONS, params.DT);
        }

        // Check boundaries:
        // =====================================================================
        //particleBC.checkBoundaryAndFlag(&params,&CS,&EB,&IONS);


        // Calculate new particle weight:
        // =====================================================================
        //particleBC.calculateParticleWeight(&params,&CS,&EB,&IONS);
        // - Count how many left the domain
        // - Calculate a_new

        // Re-inject particles:
        // =====================================================================
        //particleBC.applyParticleReinjection(&params,&CS,&EB,&IONS);
        // Loop over all Particles
        // Use f1 and f2 flags and assigne a = a_new
        // RE-inject particle states based on plasma source type
        // reset f1 and f2

        // Calculate ion moments:
        // =====================================================================
        //ionsDynamics.extrapolateIonsMoments(&params, &EB, &IONS);
        // - Apply the "a" on the extrapolation but not interpolation.


        // Apply collision operator:
        // =====================================================================
        if (params.SW.Collisions == 1)
        {
            //FPCOLL.ApplyCollisionOperator(&params,&CS,&IONS);
        }

        // Field solve:
        // =====================================================================
        // Magnetic field:
        if (params.SW.BfieldSolve == 1)
        {
            // Use Faraday's law to advance the magnetic field to level B^(N+1).
            //fields_solver.advanceBField(&params, &EB, &IONS);
        }

        if (params.SW.EfieldSolve == 1)
        {
            // Electric field:
            if(tt > 2)
            {
                // We use the generalized Ohm's law to advance in time the Electric field to level E^(N+1).
                // Using the Bashford-Adams extrapolation.
                //fields_solver.advanceEField(&params, &EB, &IONS, true, true);
            }
            else
            {
                // Using basic velocity extrapolation:
                //fields_solver.advanceEField(&params, &EB, &IONS, true, false);
            }
        }

        // Advance time:
        // =====================================================================
        //currentTime += params.DT*CS.time;


        // Save data:
        // =====================================================================
        /*
        if(fmod((double)(tt + 1), params.outputCadenceIterations) == 0)
        {
            vector<IT> IONS_OUT = IONS;

            // The ions' velocity is advanced in time in order to obtain V^(N+1):
            ionsDynamics.advanceIonsVelocity(&params, &CS, &EB, &IONS_OUT, 0.5*params.DT);

            hdfObj.saveOutputs(&params, &IONS_OUT, &EB, &CS, outputIterator+1, currentTime);

            outputIterator++;
        }
        */

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

    cout << "Hello from PICOS++!!" << endl;

    // Finalizing MPI communications:
    // =========================================================================
    mpi_main.finalizeCommunications(&params);

    cout << "end" << endl;

	return(0);
}
