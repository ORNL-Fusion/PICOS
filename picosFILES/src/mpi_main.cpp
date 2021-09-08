#include "mpi_main.h"

void MPI_MAIN_TYP::createMPITopology(params_TYP * params)
{
    // Generation of communicators for fields and particles:
    // =====================================================
    // MPI_DOMAIN_NUMBER is defined in the constructor of init_TYP

	// Define the color of each MPI sub-communicator:
	if (params->mpi.MPI_DOMAIN_NUMBER < params->mpi.MPIS_FIELDS)
    {
        params->mpi.COMM_COLOR = FIELDS_MPI_COLOR; // Color of fields' communicator
    }
    else
    {
        params->mpi.COMM_COLOR = PARTICLES_MPI_COLOR; // Color of particles' communicator
    }

	// Generate communicators for fields and particle MPI processes:
	// key = params->mpi.MPI_DOMAIN_NUMBER
    MPI_Comm_split(MPI_COMM_WORLD, params->mpi.COMM_COLOR, params->mpi.MPI_DOMAIN_NUMBER, &params->mpi.COMM);

	// Get RANK and SIZE in the context of the new communicators:
    MPI_Comm_rank(params->mpi.COMM, &params->mpi.COMM_RANK);
	MPI_Comm_size(params->mpi.COMM, &params->mpi.COMM_SIZE);

	// Define PARTICLES master nodes and distribute info among all MPIs in MPI_COMM_WORLD:
	// ==================================================================================
	// Indentify PARTICLES master node:
    if ((params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR) && (params->mpi.COMM_RANK == 0))
    {
		params->mpi.IS_PARTICLES_ROOT = true;
		params->mpi.PARTICLES_ROOT_WORLD_RANK = params->mpi.MPI_DOMAIN_NUMBER;
	}
    else
    {
		params->mpi.IS_PARTICLES_ROOT = false;
		params->mpi.PARTICLES_ROOT_WORLD_RANK = -1;
	}

	// Broadcast PARTICLES_ROOT_WORLD_RANK fom PARTICLE root to all processes in PARTICLE COMM:
	if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
	{
		MPI_Bcast(&params->mpi.PARTICLES_ROOT_WORLD_RANK, 1, MPI_INT, 0, params->mpi.COMM);
	}

	// Broadcast from a PARTICLE node to ALL nodes in COMM_WORLD:
	MPI_Bcast(&params->mpi.PARTICLES_ROOT_WORLD_RANK, 1, MPI_INT, (params->mpi.NUMBER_MPI_DOMAINS - 1), MPI_COMM_WORLD);

	// Define FIELDS master nodes and distribute info among all MPIs in MPI_COMM_WORLD:
	// ================================================================================
	// Indentify FIELDS master node:
	if ((params->mpi.COMM_COLOR == FIELDS_MPI_COLOR) && (params->mpi.COMM_RANK == 0))
    {
		params->mpi.IS_FIELDS_ROOT = true;
		params->mpi.FIELDS_ROOT_WORLD_RANK = params->mpi.MPI_DOMAIN_NUMBER;
	}
    else
    {
		params->mpi.IS_FIELDS_ROOT = false;
		params->mpi.FIELDS_ROOT_WORLD_RANK = -1;
	}

	// Broadcast FIELDS_ROOT_WORLD_RANK fom FIELDS root to all processes in FIELDS COMM:
	if (params->mpi.COMM_COLOR == FIELDS_MPI_COLOR)
    {
		MPI_Bcast(&params->mpi.FIELDS_ROOT_WORLD_RANK, 1, MPI_INT, 0, params->mpi.COMM);
    }

	// Broadcast from a FIELDS node to ALL nodes in COMM_WORLD:
	MPI_Bcast(&params->mpi.FIELDS_ROOT_WORLD_RANK, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// cout << "WORLD RANK: " << params->mpi.MPI_DOMAIN_NUMBER << " | IS PM: " << params->mpi.IS_PARTICLES_ROOT << " | IS FM: " << params->mpi.IS_FIELDS_ROOT << " | PMWR: " <<  params->mpi.PARTICLES_ROOT_WORLD_RANK << " | FMWR:" << params->mpi.FIELDS_ROOT_WORLD_RANK << endl;

    // Cartesian topology generation for MPIs evolving the fields in the simulation:
	// ============================================================================
    if (params->mpi.COMM_COLOR == FIELDS_MPI_COLOR)
    {
        int ndims;
        int dims_1D[1] = {params->mpi.MPIS_FIELDS};
        int reorder(0);
        int periods_1D[1] = {1};
        int src;
        int coord;
        int coords_1D[1];
        int topo_status;

        ndims = 1;
        params->mpi.MPI_DOMAINS_ALONG_X_AXIS = params->mpi.MPIS_FIELDS;

        MPI_Cart_create(params->mpi.COMM, ndims, dims_1D, periods_1D, reorder, &params->mpi.MPI_TOPO);

        MPI_Topo_test(params->mpi.MPI_TOPO, &topo_status);

        if (topo_status == MPI_CART)
        {
            MPI_Comm_rank(params->mpi.MPI_TOPO, &params->mpi.MPI_DOMAIN_NUMBER_CART);

            MPI_Cart_coords(params->mpi.MPI_TOPO, params->mpi.MPI_DOMAIN_NUMBER_CART, ndims, params->mpi.MPI_CART_COORDS_1D);

			for (int mpis=0; mpis<params->mpi.MPIS_FIELDS; mpis++)
			{
				int * COORDS = new int[2*params->mpi.MPIS_FIELDS];

				MPI_Allgather(&params->mpi.MPI_CART_COORDS_2D, 2, MPI_INT, COORDS, 2, MPI_INT, params->mpi.MPI_TOPO);

				params->mpi.MPI_CART_COORDS.push_back(new int[2]);
				*(params->mpi.MPI_CART_COORDS.at(mpis)) = *(COORDS + 2*mpis);
				*(params->mpi.MPI_CART_COORDS.at(mpis) + 1) = *(COORDS + 2*mpis + 1);

				delete[] COORDS;
			}

			// Calculate neighboring RANKS:
			// ===========================
			// Right side:
			src = params->mpi.MPI_DOMAIN_NUMBER_CART;
			MPI_Cart_shift(params->mpi.MPI_TOPO, 0, 1, &src, &params->mpi.RIGHT_MPI_DOMAIN_NUMBER_CART);

			// Left side:
			src = params->mpi.MPI_DOMAIN_NUMBER_CART;
			MPI_Cart_shift(params->mpi.MPI_TOPO, 0, -1, &src, &params->mpi.LEFT_MPI_DOMAIN_NUMBER_CART);

			// Calculate mesh index range (iIndex,fIndex) associated with each RANK:
            params->mpi.iIndex = params->mesh.NX_PER_MPI*params->mpi.MPI_DOMAIN_NUMBER_CART+1;
            params->mpi.fIndex = params->mesh.NX_PER_MPI*(params->mpi.MPI_DOMAIN_NUMBER_CART+1);

			/*
			cout << "World Rank: " << params->mpi.MPI_DOMAIN_NUMBER <<" , Cart Rank: "<< params->mpi.MPI_DOMAIN_NUMBER_CART << endl;
			cout << "World Rank: " << params->mpi.MPI_DOMAIN_NUMBER <<" , iIndex: "<< params->mpi.iIndex << endl;
			cout << "World Rank: " << params->mpi.MPI_DOMAIN_NUMBER <<" , fIndex: "<< params->mpi.fIndex << endl;

			World Rank: 0 , Cart Rank: 0
			World Rank: 0 , iIndex: 1
			World Rank: 0 , fIndex: 250
			World Rank: 1 , Cart Rank: 1
			World Rank: 1 , iIndex: 251
			World Rank: 1 , fIndex: 500
			*/

			if (params->mpi.MPI_DOMAIN_NUMBER_CART == 0)
            {
                cout << endl << "* * * * * * * * * * * * GENERATING MPI TOPOLOGY FOR FIELDS * * * * * * * * * * * * * * * * * *" << endl;
                cout << "+ Number of MPI processes along the x-axis: " << dims_1D[0] << endl;
                cout << "+ Number of mesh nodes along x-axis: " << params->mesh.NX_IN_SIM << endl;
				cout << "+ Number of MPI processes for FIELDS: " << params->mpi.MPIS_FIELDS << endl;
				cout << "+ Number of MPI processes for PARTICLES: " << params->mpi.MPIS_PARTICLES << endl;
				cout << "* * * * * * * * * * * * MPI TOPOLOGY FOR FIELDS GENERATED  * * * * * * * * * * * * * * * * * *" << endl << endl;
			}
		}
        else
        {
			cerr << "ERROR: MPI topology could not be created!" << endl;

			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Abort(MPI_COMM_WORLD,-102);
		}
	}

}

// Implement methods: finalizeCommunications
// ==================================================================================================================================
void MPI_MAIN_TYP::finalizeCommunications(params_TYP * params)
{
	bool finalized = false;

	if(params->mpi.MPI_DOMAIN_NUMBER == 0)
		cout << endl << "* * * * * * * * * * * * FINALIZING MPI COMMUNICATIONS * * * * * * * * * * * * * * * * * *" << endl;

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Comm_free(&params->mpi.COMM);

	MPI_Finalize();

	int temp;
	MPI_Finalized(&temp);
	finalized = (bool)temp;

	if(finalized)
		cout << "MPI process: " << params->mpi.MPI_DOMAIN_NUMBER << " FINALIZED" << endl;
	else
		cout << "MPI process: " << params->mpi.MPI_DOMAIN_NUMBER << " NOT FINALIZED - ERROR" << endl;

	if(params->mpi.MPI_DOMAIN_NUMBER == 0)
		cout << "* * * * * * * * * * * * MPI COMMUNICATIONS FINALIZED * * * * * * * * * * * * * * * * * *" << endl;
}
