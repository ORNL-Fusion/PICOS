#include "fieldSolve.h"

// Constructor:
// ============================================================================================
fields_solver_TYP::fields_solver_TYP(const params_TYP * params, CS_TYP * CS)
{
	NX_S = params->mesh.NX_PER_MPI + 2; // Grid cells per domain + 2 ghost cells
	NX_T = params->mesh.NX_IN_SIM + 2;  // Grid cells in entire simulation + 2 ghost cells
	NX_R = params->mesh.NX_IN_SIM;      // Grid cells in entire simulation

    // Electron density at various time steps:
	ne.zeros(NX_S);
	ne_.zeros(NX_S);
	ne__.zeros(NX_S);
	ne___.zeros(NX_S);

	// Electron temperature:
	Te.zeros(NX_S);

	// Electron pressure:
	Pe.zeros(NX_S);

	// Electron pressure gradient:
	dPe.zeros(NX_S);

	//nU.zeros(NX_S);
    //U.zeros(NX_S);

    // Electron density gradient:
	//dne.zeros(NX_S);
	dx = params->mesh.DX;

    // Electric field:
    EX_m.zeros(NX_S);
}

// Fill ghost cells:
// ============================================================================================
void fields_solver_TYP::fillGhosts(arma::vec * C)
{
	int NX = C->n_elem;

	(*C)(0)    = (*C)(1);
	(*C)(NX-1) = (*C)(NX-2);
}

void fields_solver_TYP::fill4Ghosts(arma::vec * v)
{
	int NX = v->n_elem;

	v->subvec(0,1)       = v->subvec(2,3);
	v->subvec(NX-2,NX-1) = v->subvec(NX-4,NX-3);
}

// Smoothing:
// ============================================================================================
void fields_solver_TYP::smooth(arma::vec * v, double as)
{
	int NX = v->n_elem;

	arma::vec b = zeros(NX);

	double wc(0.5); 	// center weight
	double ws(0.25);	// sides weight

	//Step 1: Averaging process
	b.subvec(1, NX-2) = v->subvec(1, NX-2);

	fillGhosts(&b);

	b.subvec(1, NX-2) = wc*b.subvec(1, NX-2) + ws*b.subvec(2, NX-1) + ws*b.subvec(0, NX-3);

	//Step 2: Averaged weighted variable estimation.
	v->subvec(1, NX-2) = (1.0 - as)*v->subvec(1, NX-2) + as*b.subvec(1, NX-2);
}

// MPI functions:
// ============================================================================================
void fields_solver_TYP::MPI_Allgathervec(const params_TYP * params, arma::vec * field)
{
    // Rank-dependent subdomain indices:
    unsigned int iIndex = params->mpi.iIndex;
    unsigned int fIndex = params->mpi.fIndex;

    // Buffers:
    arma::vec recvbuf(params->mesh.NX_IN_SIM);
    arma::vec sendbuf(params->mesh.NX_PER_MPI);

    // Allgather for x-component
    sendbuf = field->subvec(iIndex, fIndex);
    MPI_Allgather(sendbuf.memptr(), params->mesh.NX_PER_MPI, MPI_DOUBLE, recvbuf.memptr(), params->mesh.NX_PER_MPI, MPI_DOUBLE, params->mpi.MPI_TOPO);

    // Assign to output:
    field->subvec(1, params->mesh.NX_IN_SIM) = recvbuf;
}

void fields_solver_TYP::MPI_SendVec(const params_TYP * params, arma::vec * v)
{
	// Send vector v from FIELDS root to PARTICLES root:
	if (params->mpi.IS_FIELDS_ROOT)
	{
		MPI_Send(v->memptr(), v->n_elem, MPI_DOUBLE, params->mpi.PARTICLES_ROOT_WORLD_RANK, FIELDS_TAG, MPI_COMM_WORLD);
	}

    // Receive vector v from FIELDS root into PARTICLES root:
	if (params->mpi.IS_PARTICLES_ROOT)
	{
		MPI_Recv(v->memptr(), v->n_elem, MPI_DOUBLE, params->mpi.FIELDS_ROOT_WORLD_RANK, FIELDS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// PARTICLES root broadcasts to all PARTICLE ranks:
	if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
	{
		MPI_Bcast(v->memptr(), v->n_elem, MPI_DOUBLE, 0, params->mpi.COMM);
	}
}

// Electric field solve:
// ============================================================================================
void fields_solver_TYP::advanceEfield(const params_TYP * params, fields_TYP * fields, CS_TYP * CS, vector<ionSpecies_TYP> * IONS, electrons_TYP * electrons)
{
	if (params->mpi.COMM_COLOR == FIELDS_MPI_COLOR)
	{
		// Indices of subdomain:
		unsigned int iIndex = params->mpi.iIndex;
		unsigned int fIndex = params->mpi.fIndex;

		// Initialize the electron density:
		ne.zeros();
		ne_.zeros();
		ne__.zeros();
		ne___.zeros();

		// Initialize the electron temperature:
		//Te.zeros();

		// Initialize the electron pressure:
		//Pe.zeros();

        // Initialize the gradient:
		//dPe.zeros();
        //dne.zeros();

		// Electron temperature:
		//double Te = params->f_IC.Te;

		// Accumulate the contribution from all species:
		for(int ss=0; ss<IONS->size(); ss++)
		{
			double Z = IONS->at(ss).Z;

			// Electron density:
			ne    += Z*IONS->at(ss).n_m.subvec(iIndex - 1, fIndex + 1);
			ne_   += Z*IONS->at(ss).n_m_.subvec(iIndex - 1, fIndex + 1);
			ne__  += Z*IONS->at(ss).n_m__.subvec(iIndex - 1, fIndex + 1);
			ne___ += Z*IONS->at(ss).n_m___.subvec(iIndex - 1, fIndex + 1);

		}

		// Time-averaged electron density;
		ne = (ne + ne_ + ne__ + ne___)/4;
        fill4Ghosts(&ne);

		// Electron temperature:
		Te = electrons->Te_m.subvec(iIndex - 1, fIndex + 1);

		// Electron pressure:
		Pe = (Te%ne)/F_E_DS;

		// Gradient in the electron pressure:
		dPe.subvec(1,NX_S - 2) = 0.5*( Pe.subvec(2,NX_S-1) - Pe.subvec(0,NX_S-3) );
        fill4Ghosts(&dPe);

		// Gradient of electron density:
		//dne.subvec(1,NX_S - 2) = 0.5*( ne.subvec(2,NX_S-1) - ne.subvec(0,NX_S-3) );
        //fill4Ghosts(&dne);

		// Electric field based on Ohm's law:
        // EX_m = -(Te/F_E_DS)*(1/ne)%dne/dx;
		EX_m = -(1/ne)%dPe/dx;
		fields->EX_m.subvec(iIndex,fIndex) = EX_m.subvec(1,NX_S - 2);
        fill4Ghosts(&fields->EX_m);

		// Error checks:
		// =============
		#ifdef CHECKS_ON
		if(!fields->EX_m.is_finite())
		{
			cout << "Non finite values in Ex" << endl;
			MPI_Abort(params->mpi.MPI_TOPO, -110);
		}
		#endif

		// Assemble entire profile:
		// =======================
		MPI_Allgathervec(params, &fields->EX_m);

		// Apply smoothing:
        // ===============
        for (int jj=0; jj<params->filtersPerIterationFields; jj++)
        {
			smooth(&fields->EX_m, params->smoothingParameter);
		}

	} // FIELDS MPI

	// Send to PARTICLE ranks:
	// ======================
	MPI_SendVec(params,&fields->EX_m);
}
