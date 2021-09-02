#include "PIC.h"

void PIC_TYP::MPI_AllreduceVec(const params_TYP * params, arma::vec * v)
{
	arma::vec recvbuf = zeros(v->n_elem);

	MPI_Allreduce(v->memptr(), recvbuf.memptr(), v->n_elem, MPI_DOUBLE, MPI_SUM, params->mpi.MPI_TOPO);

	*v = recvbuf;
}

void PIC_TYP::MPI_SendVec(const params_TYP * params, arma::vec * v)
{
	// Send vector from PARTICLE ROOT to FIELDS ROOT:
	if (params->mpi.IS_PARTICLES_ROOT)
    {
		MPI_Send(v->memptr(), v->n_elem, MPI_DOUBLE, params->mpi.FIELDS_ROOT_WORLD_RANK, PARTICLES_TAG, MPI_COMM_WORLD);
	}

	if (params->mpi.IS_FIELDS_ROOT)
    {
		MPI_Recv(v->memptr(), v->n_elem, MPI_DOUBLE, params->mpi.PARTICLES_ROOT_WORLD_RANK, PARTICLES_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// Broadcast from FIELDS ROOT to all other FIELDS MPIs:
	if (params->mpi.COMM_COLOR == FIELDS_MPI_COLOR)
    {
        MPI_Bcast(v->memptr(), v->n_elem, MPI_DOUBLE, 0, params->mpi.COMM);
    }
}


void PIC_TYP::MPI_ReduceVec(const params_TYP * params, arma::vec * v)
{
    // Create receive buffer:
    // ======================
    arma::vec recvbuf = zeros(v->n_elem);

    // Reduce the vector at the PARTICLE ROOT:
    // ======================================
    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        MPI_Reduce(v->memptr(), recvbuf.memptr(), v->n_elem, MPI_DOUBLE, MPI_SUM, 0, params->mpi.COMM);

        if (params->mpi.IS_PARTICLES_ROOT)
        {
             *v = recvbuf;
        }
    }
}

void PIC_TYP::MPI_Allgathervec(const params_TYP * params, arma::vec * field)
{
	unsigned int iIndex = params->mpi.iIndex;
	unsigned int fIndex = params->mpi.fIndex;

	arma::vec recvbuf(params->mesh.NX_IN_SIM);
	arma::vec sendbuf(params->mesh.NX_PER_MPI);

	sendbuf = field->subvec(iIndex, fIndex);
	MPI_Allgather(sendbuf.memptr(), params->mesh.NX_PER_MPI, MPI_DOUBLE, recvbuf.memptr(), params->mesh.NX_PER_MPI, MPI_DOUBLE, params->mpi.MPI_TOPO);
	field->subvec(1, params->mesh.NX_IN_SIM) = recvbuf;
}


void PIC_TYP::MPI_Recvvec(const params_TYP * params, arma::vec * field)
{
	// We send the vector from root process of fields to root process of particles
	arma::vec recvbuf(params->mesh.NX_IN_SIM);
	arma::vec sendbuf(params->mesh.NX_IN_SIM);

	sendbuf = field->subvec(1, params->mesh.NX_IN_SIM);

	if (params->mpi.IS_FIELDS_ROOT)
    {
		MPI_Send(sendbuf.memptr(), params->mesh.NX_IN_SIM, MPI_DOUBLE, params->mpi.PARTICLES_ROOT_WORLD_RANK, 0, MPI_COMM_WORLD);
	}

	if (params->mpi.IS_PARTICLES_ROOT)
    {
		MPI_Recv(recvbuf.memptr(), params->mesh.NX_IN_SIM, MPI_DOUBLE, params->mpi.FIELDS_ROOT_WORLD_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		field->subvec(1, params->mesh.NX_IN_SIM) = recvbuf;
	}

	// Then, the fields is broadcasted to all processes in the particles communicator COMM
	if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
		sendbuf = field->subvec(1, params->mesh.NX_IN_SIM);

		MPI_Bcast(sendbuf.memptr(), params->mesh.NX_IN_SIM, MPI_DOUBLE, 0, params->mpi.COMM);

		field->subvec(1, params->mesh.NX_IN_SIM) = sendbuf;
	}
}

void PIC_TYP::MPI_Recv_AllFields(const params_TYP * params, fields_TYP * fields)
{
	// Send field data from FIELDS ranks and recieve fields data at PARTICLE ranks
	MPI_Recvvec(params,&fields->EX_m);
	MPI_Recvvec(params,&fields->BX_m);
	MPI_Recvvec(params,&fields->dBX_m);

	if (params->SW.RFheating == 1)
	{
		MPI_Recvvec(params,&fields->ddBX_m);
	}
}

void PIC_TYP::fillGhosts(arma::vec * C)
{
	int NX = C->n_elem;

	(*C)(0)    = (*C)(1);
	(*C)(NX-1) = (*C)(NX-2);
}

void PIC_TYP::fillGhost_AllFields(const params_TYP * params, fields_TYP * fields)
{
	fillGhosts(&fields->EX_m);
	fillGhosts(&fields->BX_m);
	fillGhosts(&fields->dBX_m);

	if (params->SW.RFheating == 1)
	{
		fillGhosts(&fields->ddBX_m);
	}
}

void PIC_TYP::fill4Ghosts(arma::vec * v)
{
	int NX = v->n_elem;

	v->subvec(0,1)       = v->subvec(2,3);
	v->subvec(NX-2,NX-1) = v->subvec(NX-4,NX-3);
}

void PIC_TYP::smooth(arma::vec * v, double as)
{
	int NX(v->n_elem);

	arma::vec b = zeros(NX);

	double wc(0.75); 	// center weight
	double ws(0.125);	// sides weight

	//Step 1: Averaging process
	b.subvec(1, NX-2) = v->subvec(1, NX-2);

	fillGhosts(&b);

	b.subvec(1, NX-2) = wc*b.subvec(1, NX-2) + ws*b.subvec(2, NX-1) + ws*b.subvec(0, NX-3);

	//Step 2: Averaged weighted variable estimation.
	v->subvec(1, NX-2) = (1.0 - as)*v->subvec(1, NX-2) + as*b.subvec(1, NX-2);
}

// Constructor:
PIC_TYP::PIC_TYP(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
	// Get latest mesh-defined values from FIELDS ranks:
	// =================================================
	MPI_Recv_AllFields(params,fields);

	// Fill the ghost cells in all fields:
	fillGhost_AllFields(params,fields);

	// Assign cell for all particles:
    // =============================
	assignCell_AllSpecies(params,IONS);

	// Interpolate all fields on all species:
	// ======================================
	interpolateFields_AllSpecies(params,IONS,fields);

	// Calculate ion moments and populate mesh-defined ion moments:
	// ============================================================
	// Run 3 dummy cycles to load "n" and "nv" at previous time steps:
	for(int tt=0; tt<3; tt++)
	{
		extrapolateMoments_AllSpecies(params,CS,fields,IONS);
	}

}

void PIC_TYP::interpolateScalarField(const params_TYP * params, ionSpecies_TYP * IONS, const arma::vec * F_m, arma::vec * F_p)
{
    int NX =  params->mesh.NX_IN_SIM + 4; //Mesh size along the X axis (considering the gosht cell)
    int NSP(IONS->NSP);

    // Allocate memory and initialize to zero:
    arma::vec F = zeros(NX);

    // Fill in vector F with mesh-defined data:
    F.subvec(1,NX-2) = *F_m;

    // Take care of ghost cells:
    fill4Ghosts(&F);

    #pragma omp parallel for default(none) shared(params, IONS, F_p, F) firstprivate(NSP)
    for(int ii=0; ii<NSP; ii++)
    {
        int ix = IONS->mn(ii) + 2;

        (*F_p)(ii) += IONS->wxl(ii)*F(ix-1);
        (*F_p)(ii) += IONS->wxc(ii)*F(ix);
        (*F_p)(ii) += IONS->wxr(ii)*F(ix+1);

    }// omp parallel for
}

void PIC_TYP::interpolateFields(const params_TYP * params, ionSpecies_TYP * IONS, const fields_TYP * fields)
{
	// Need to use SW in order to enable/disable ddBX interpolation

	// Allocate memory for field variables:
	arma::vec EX_p   = zeros(IONS->NSP, 1);
	arma::vec BX_p   = zeros(IONS->NSP, 1);
	arma::vec dBX_p  = zeros(IONS->NSP, 1);

	// Interpolate mesh-defined fields into particles:
	interpolateScalarField(params, IONS, &fields->EX_m , &EX_p );
	interpolateScalarField(params, IONS, &fields->BX_m , &BX_p );
	interpolateScalarField(params, IONS, &fields->dBX_m, &dBX_p);

	// Assign values:
	IONS->EX_p   = EX_p;
	IONS->BX_p   = BX_p;
	IONS->dBX_p  = dBX_p;

	if (params->SW.RFheating == 1)
	{
		arma::vec ddBX_p = zeros(IONS->NSP, 1);
		interpolateScalarField(params, IONS, &fields->ddBX_m, &ddBX_p);
		IONS->ddBX_p = ddBX_p;
	}
}

void PIC_TYP::interpolateFields_AllSpecies(const params_TYP * params, vector<ionSpecies_TYP> * IONS, const fields_TYP * fields)
{
	// Interpolate all fields on all species:
	// =======================
	for(int ss=0;ss<IONS->size();ss++)
	{
		if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
		{
			// Interpolate mesh-defined fields into ALL particles locations:
			interpolateFields(params, &IONS->at(ss), fields);
		}
	}
}

void PIC_TYP::interpEM(const params_TYP * params, const ionSpecies_TYP * IONS, const fields_TYP * fields, arma::rowvec * ZN, arma::rowvec * EM)
{
    // Assign cell:
    double xp   = (*ZN)(0);
    double xMin = 0;
    double DX   = params->mesh.DX;
    int m       = round( 0.5 + (xp - xMin)/DX ) - 1;

    // During RK4, some projections can be go out of bound
    if ( m >= params->mesh.NX_IN_SIM)
    {
        m = params->mesh.NX_IN_SIM - 1;
    }
    if ( m < 0)
    {
        m = 0;
    }

    // Distance to nearest grid point:
    double X    = params->mesh.nodesX(m) - xp;

    // Assignment function:
    arma::vec W = zeros<vec>(3);
    W(0) = 0.5*pow( 1.5 + ((X - DX)/DX) ,2); // Left:
    W(1) = 0.75 - pow(X/DX,2);               // Center:
    W(2) = 0.5*pow(1.5 - ((X + DX)/DX) ,2 ); // Right:

    // Nearest grid point:
    int ix = m + 1;

    // Temporary storage for interpolated fields:
    arma::vec f = zeros<vec>(3);

    // Interpolate:
    // EX:
    f(0) = fields->EX_m(ix - 1);
    f(1) = fields->EX_m(ix);
    f(2) = fields->EX_m(ix + 1);
    (*EM)(0) = arma::dot(f,W);

    // BX:
    f(0) = fields->BX_m(ix - 1);
    f(1) = fields->BX_m(ix);
    f(2) = fields->BX_m(ix + 1);
    (*EM)(1) = arma::dot(f,W);

    // dBX:
    f(0) = fields->dBX_m(ix - 1);
    f(1) = fields->dBX_m(ix);
    f(2) = fields->dBX_m(ix + 1);
    (*EM)(2) = arma::dot(f,W);

}

void PIC_TYP::calculateF(const params_TYP * params, const ionSpecies_TYP * IONS, arma::rowvec * ZN, arma::rowvec * EM, arma::rowvec * F)
{
    // Ion parameters:
    double qa = IONS->Q;
    double Ma = IONS->M;

	// Gather fields:
    double E    = (*EM)(0);
    double B    = (*EM)(1);
    double dB   = (*EM)(2);

    // Output:
	if (params->advanceParticleMethod == 1)
	{
		// Gather particle states:
		double vpar = (*ZN)(1);
		double vper = (*ZN)(2);

		// Output:
	    (*F)(0) = +vpar;
		(*F)(1) = -0.5*vper*vper*dB/B + (qa/Ma)*E;
    	(*F)(2) = +0.5*vper*vpar*dB/B;
	}
	if (params->advanceParticleMethod == 2)
	{
		// Gather particle states:
		double vpar = (*ZN)(1);
		double mu = (*ZN)(2);

	    (*F)(0) = +vpar;
		(*F)(1) = -(mu/Ma)*dB + (qa/Ma)*E;
    	(*F)(2) = 0;
	}
}

void PIC_TYP::advanceParticles(const params_TYP * params, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    // Get latest mesh-defined values from FIELDS MPIs:
	MPI_Recv_AllFields(params,fields);

	// Fill the ghost cells in all fields:
	fillGhost_AllFields(params,fields);

	// Iterate over all the ion species:
	for(int ss=0; ss<IONS->size(); ss++)
	{
		if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
		{
			// Number of particles:
			int NSP = IONS->at(ss).NSP;

            // Time step:
            double DT = params->DT;

			// Ion mass:
			double Ma = IONS->at(ss).M;

			#pragma omp parallel for default(none) shared(IONS, params, fields, DT, ss, Ma, std::cout) firstprivate(NSP, F_C_DS)
            for(int ii=0;ii<NSP;ii++)
			{
                // Start RK4 solution:
                //==============================================================
                /*
                ! Equations of motion are expressed as follows:
                !
                ! dZ/dt = F, thus Z1 = Z0 + dZ
                !
                ! where
                !
                ! dZ = F(Z0)*dt
                ! Z = [z, vz, vp]
                ! F(1) = +vz
                ! F(2) = -0.5*vp*vp*dB/B + (q/Ma)*E
                ! F(3) = +0.5*vp*vz*dB/B
                */

				// Assemble vectors to use in RK4 method:
				arma::rowvec Z0  = zeros<rowvec>(3);
                arma::rowvec ZN  = zeros<rowvec>(3);
				arma::rowvec Z1  = zeros<rowvec>(3);
				arma::rowvec EM  = zeros<rowvec>(3);
                arma::rowvec F   = zeros<rowvec>(3);
                arma::rowvec dZ1 = zeros<rowvec>(3);
                arma::rowvec dZ2 = zeros<rowvec>(3);
                arma::rowvec dZ3 = zeros<rowvec>(3);
                arma::rowvec dZ4 = zeros<rowvec>(3);

				// Extract particle states:
				double x    = IONS->at(ss).X_p(ii);
				double vpar = IONS->at(ss).V_p(ii,0);
				double vper = IONS->at(ss).V_p(ii,1);

				// Extract particle-defined fields:
				double E  = IONS->at(ss).EX_p(ii);
				double B  = IONS->at(ss).BX_p(ii);
				double dB = IONS->at(ss).dBX_p(ii);

				// Initialize initial particle state Z0:
				Z0 = {x, vpar, vper};

				// Interpolate fields at current particle postion:
				interpEM(params, &IONS->at(ss), fields, &Z0, &EM);

				// Select solution method:
				switch (params->advanceParticleMethod)
				{
				case 1:
					// Do nothing
					break;
				case 2:
					// Use magnetic moment
					B  = EM(1);
					double mu = 0.5*Ma*pow(vper,2)/B;
					Z0(2) = mu;
					break;
				}

                // Step 1:
                ZN = Z0;
                calculateF(params, &IONS->at(ss), &ZN, &EM, &F);
                dZ1 = F*DT;

				if ( isnan(ZN(0)) || isnan(ZN(1)) || isnan(ZN(2)) )
				{
					cout << "Z0(0), 1 = " << Z0(0) << endl;
					cout << "ZN(0), 1 = " << ZN(0) << endl;
					cout << "ZN(1), 1 = " << ZN(1) << endl;
					cout << "ZN(2), 1 = " << ZN(2) << endl;

					cout << "E  = " << E  << endl;
					cout << "B  = " << B  << endl;
					cout << "dB = " << dB << endl;
					cout << "EM(0) = " << EM(0) << endl;
					cout << "EM(1) = " << EM(1) << endl;
					cout << "EM(2) = " << EM(2) << endl;
				}

                // Step 2:
                ZN = Z0 + dZ1/2;
                interpEM(params, &IONS->at(ss), fields, &ZN, &EM);
                calculateF(params, &IONS->at(ss), &ZN, &EM, &F);
                dZ2 = F*DT;

				if ( isnan(ZN(0)) || isnan(ZN(1)) || isnan(ZN(2)) )
				{
					cout << "Z0(0), 2 = " << Z0(0) << endl;
					cout << "ZN(0), 2 = " << ZN(0) << endl;
					cout << "ZN(1), 2 = " << ZN(1) << endl;
					cout << "ZN(2), 2 = " << ZN(2) << endl;

					cout << "EM(0), 2 = " << EM(0) << endl;
					cout << "EM(1), 2 = " << EM(1) << endl;
					cout << "EM(2), 2 = " << EM(2) << endl;
				}

                // Step 3:
                ZN = Z0 + dZ2/2;
                interpEM(params, &IONS->at(ss), fields, &ZN, &EM);
                calculateF(params, &IONS->at(ss), &ZN, &EM, &F);
                dZ3 = F*DT;

				if ( isnan(ZN(0)) || isnan(ZN(1)) || isnan(ZN(2)) )
				{
					cout << "Z0(0), 3 = " << Z0(0) << endl;
					cout << "ZN(0), 3 = " << ZN(0) << endl;
					cout << "ZN(1), 3 = " << ZN(1) << endl;
					cout << "ZN(2), 3 = " << ZN(2) << endl;
				}

                // Step 4:
                ZN = Z0 + dZ3;
                interpEM(params, &IONS->at(ss), fields, &ZN, &EM);
                calculateF(params, &IONS->at(ss), &ZN, &EM, &F);
                dZ4 = F*DT;

				if ( isnan(ZN(0)) || isnan(ZN(1)) || isnan(ZN(2)) )
				{
					cout << "Z0(0), 4 = " << Z0(0) << endl;
					cout << "ZN(0), 4 = " << ZN(0) << endl;
					cout << "ZN(1), 4 = " << ZN(1) << endl;
					cout << "ZN(2), 4 = " << ZN(2) << endl;
				}

                // Assemble RK4 solution:
                Z1 = Z0 + (dZ1 + 2*dZ2 + 2*dZ3 + dZ4)/6;

				if ( isnan(ZN(0)) || isnan(ZN(1)) || isnan(ZN(2)) )
				{
					cout << "Z0(0) = " << Z0(0) << endl;
					cout << "Z1(0) = " << Z1(0) << endl;
					cout << "Z1(1) = " << Z1(1) << endl;
					cout << "Z1(2) = " << Z1(2) << endl;
				}

				// Interpolate fields at new particle position:
				interpEM(params, &IONS->at(ss), fields, &Z1, &EM);

				// Assign solution to output vector:
				switch (params->advanceParticleMethod)
				{
				case 1:
					// Do nothing, RK4 solved for (x, vpar, vper)
					break;
				case 2:
					// Calculate vper since RK4 solved for (x, vpar, mu)
					B  = EM(1);
					double mu = Z1(2);
					double vper = sqrt(2*mu*B/Ma);
					Z1(2)  = vper;
					break;
				}
                // End of RK solution:
                //==============================================================

                // Update new particle states:
                IONS->at(ss).X_p(ii)   = Z1(0);
                IONS->at(ss).V_p(ii,0) = Z1(1); // vpar
                IONS->at(ss).V_p(ii,1) = Z1(2); // vper
				IONS->at(ss).mu_p(ii)  = 0.5*Ma*pow(Z1(2),2)/EM(1) ; // mu

			} // End of parallel region
		}
	}//structure to iterate over all the ion species.
}

void PIC_TYP::assignCell(const params_TYP * params, ionSpecies_TYP * IONS)
{
	// Total number of computational particles:
	int NSP(IONS->NSP);

	// Clear assignment function:
    IONS->wxc.zeros();
    IONS->wxl.zeros();
    IONS->wxr.zeros();

	#pragma omp parallel for default(none) shared(IONS, params, std::cout) firstprivate(NSP)
    for(int ii=0; ii<NSP; ii++)
    {
		// Calculate nearest grid point:
		double X_p     = IONS->X_p(ii);
		double X_p_min = 0;
		double DX      = params->mesh.DX;
		int m = round( 0.5 + (X_p - X_p_min)/DX ) - 1;

		// Correct "m" near boundaries if out of bound:
		if ( m >= params->mesh.NX_IN_SIM)
		{
			m = params->mesh.NX_IN_SIM - 1;
			//cout << "m exceeds boundary" << endl;
			//cout << "X_p = " << X_p << endl;
			//cout << "LX = " << params->mesh.LX << endl;
		}
		if ( m < 0)
		{
			//cout << "m = " << m << endl;
			m = 0;
			//cout << "m exceeds boundary" << endl;
			//cout << "X_p = " << X_p << endl;
		}

		// Assign nearest grid point:
		IONS->mn(ii) = m;

		// Distance to nearest grid point:
		double X = params->mesh.nodesX(m) - X_p;

		// Assignment function:
		IONS->wxl(ii) = 0.5*pow(1.5 + ((X - DX)/DX),2); // Left:
		IONS->wxc(ii) = 0.75 - pow(X/DX,2);               // Center:
		IONS->wxr(ii) = 0.5*pow(1.5 - ((X + DX)/DX),2); // Right:

	} // parallel omp

}

void PIC_TYP::assignCell_AllSpecies(const params_TYP * params, vector<ionSpecies_TYP> * IONS)
{
	// Iterate over all ion species:
    // =============================
    for(int ss=0;ss<IONS->size();ss++)
    {
        // Assign cell and calculate partial ion moments:
        if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
        {
            // Assign cell:
            PIC_TYP::assignCell(params, &IONS->at(ss));
        }
	}
}

void PIC_TYP::extrapolateMoments_AllSpecies(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
	// Iterate over all ion species:
    // =============================
    for(int ss=0;ss<IONS->size();ss++)
    {

        // Assign cell and calculate partial ion moments:
        if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
        {
			//Calculate partial moments:
			calculateIonMoments(params,CS,fields,&IONS->at(ss));

			// Reduce IONS moments to PARTICLE ROOT:
			// =====================================
			MPI_ReduceVec(params, &IONS->at(ss).n_m);
			MPI_ReduceVec(params, &IONS->at(ss).nv_m);
			MPI_ReduceVec(params, &IONS->at(ss).P11_m);
			MPI_ReduceVec(params, &IONS->at(ss).P22_m);

			// Broadcast to all PARTICLE ranks:
			// ================================
			MPI_Bcast(IONS->at(ss).n_m.memptr()  , IONS->at(ss).n_m.size()  , MPI_DOUBLE, 0, params->mpi.COMM);
			MPI_Bcast(IONS->at(ss).nv_m.memptr() , IONS->at(ss).nv_m.size() , MPI_DOUBLE, 0, params->mpi.COMM);
			MPI_Bcast(IONS->at(ss).P11_m.memptr(), IONS->at(ss).P11_m.size(), MPI_DOUBLE, 0, params->mpi.COMM);
			MPI_Bcast(IONS->at(ss).P22_m.memptr(), IONS->at(ss).P22_m.size(), MPI_DOUBLE, 0, params->mpi.COMM);

			// Apply smoothing:
			// ===============
			for (int jj=0; jj<params->filtersPerIterationIons; jj++)
			{
			  smooth(&IONS->at(ss).n_m  , params->smoothingParameter);
			  smooth(&IONS->at(ss).nv_m , params->smoothingParameter);
			  smooth(&IONS->at(ss).P11_m, params->smoothingParameter);
			  smooth(&IONS->at(ss).P22_m, params->smoothingParameter);
			}

			// Calculate derived ion moments: Tpar_m, Tper_m:
			// ==============================================
			calculateDerivedIonMoments(params, CS, &IONS->at(ss));
        }

		// 0th and 1st moments at various time levels are sent to fields processes:
        // =============================================================
		// Ion density:
        MPI_SendVec(params, &IONS->at(ss).n_m);
        MPI_SendVec(params, &IONS->at(ss).n_m_);
        MPI_SendVec(params, &IONS->at(ss).n_m__);
        MPI_SendVec(params, &IONS->at(ss).n_m___);

		// For ES solution, we dont need to send these to the field ranks:
		// Ion parallel flux:
		//MPI_SendVec(params, &IONS->at(ss).nv_m);
		//MPI_SendVec(params, &IONS->at(ss).nv_m_);
		//MPI_SendVec(params, &IONS->at(ss).nv_m__);

	}
}

void PIC_TYP::calculateIonMoments(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, ionSpecies_TYP * IONS)
{
	// Ion density:
	IONS->n_m___ = IONS->n_m__;
	IONS->n_m__  = IONS->n_m_;
	IONS->n_m_   = IONS->n_m;

	// Ion flux:
	IONS->nv_m__ = IONS->nv_m_;
	IONS->nv_m_  = IONS->nv_m;

	// Calculate ion moments:
	eim(params,CS,fields,IONS);
}

void PIC_TYP::eim(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, ionSpecies_TYP * IONS)
{
	// Number of particles:
	int NSP(IONS->NSP);

	// Ion mass:
	double Ma = IONS->M;

	// Reference magnetic field:
	double B0 = params->em_IC.BX; // Maybe use the current value at the reference location

	// Clearing content of ion moments:
	// ===============================
	IONS->n_m.zeros();
	IONS->nv_m.zeros();
	IONS->P11_m.zeros();
	IONS->P22_m.zeros();

	#pragma omp parallel default(none) shared(params, IONS, B0, Ma) firstprivate(NSP)
	{
		// Create private moments:
		// ======================
		arma::vec n   = zeros(params->mesh.NX_IN_SIM + 4);
		arma::vec nv  = zeros(params->mesh.NX_IN_SIM + 4);
		arma::vec P11 = zeros(params->mesh.NX_IN_SIM + 4);
		arma::vec P22 = zeros(params->mesh.NX_IN_SIM + 4);

		// Assemble moments:
		// =================
		#pragma omp for
		for(int ii=0; ii<NSP; ii++)
		{
			// Nearest grid point:
			int ix = IONS->mn(ii) + 1;

			// Particle velocity:
			double vpar = IONS->V_p(ii,0);
			double vper = IONS->V_p(ii,1);

			// vx component:
			arma::vec phi = 2*M_PI*randu<vec>(1);
			double vy = vper*cos(phi(0));

			// Particle-defined magnetic field:
			double B = IONS->BX_p(ii);

			// Compression factor:
			double c = B/B0;

			// Particle weight:
			if (params->currentTime == 0)
			{
				IONS->a_p(ii) = 1/c;
			}
			double a = IONS->a_p(ii);

			// Density:
			n(ix-1) += IONS->wxl(ii)*a*c;
			n(ix)   += IONS->wxc(ii)*a*c;
			n(ix+1) += IONS->wxr(ii)*a*c;

			// Particle flux density:
			nv(ix-1) += IONS->wxl(ii)*a*vpar*c;
			nv(ix) 	 += IONS->wxc(ii)*a*vpar*c;
			nv(ix+1) += IONS->wxr(ii)*a*vpar*c;

			// Unitary particle flux may be needed here

			// Stress tensor P11:
			P11(ix-1) += IONS->wxl(ii)*a*Ma*pow(vpar,2)*c;
			P11(ix)   += IONS->wxc(ii)*a*Ma*pow(vpar,2)*c;
			P11(ix+1) += IONS->wxr(ii)*a*Ma*pow(vpar,2)*c;

			// Stress tensor P22:
			P22(ix-1) += IONS->wxl(ii)*a*Ma*pow(vy,2)*c;
			P22(ix)   += IONS->wxc(ii)*a*Ma*pow(vy,2)*c;
			P22(ix+1) += IONS->wxr(ii)*a*Ma*pow(vy,2)*c;
		}

		// Reduce partial moments from each thread:
		// ========================================
		#pragma omp critical (update_ion_moments)
		{
			IONS->n_m.subvec(1,params->mesh.NX_IN_SIM)   += n.subvec(2,params->mesh.NX_IN_SIM + 1);
			IONS->nv_m.subvec(1,params->mesh.NX_IN_SIM)  += nv.subvec(2,params->mesh.NX_IN_SIM + 1);
			IONS->P11_m.subvec(1,params->mesh.NX_IN_SIM) += P11.subvec(2,params->mesh.NX_IN_SIM + 1);
			IONS->P22_m.subvec(1,params->mesh.NX_IN_SIM) += P22.subvec(2,params->mesh.NX_IN_SIM + 1);
		}

	}//End of the parallel region

	// Ghost contributions:
	// ====================
	fillGhosts(&IONS->n_m);
	fillGhosts(&IONS->nv_m);
	fillGhosts(&IONS->P11_m);
	fillGhosts(&IONS->P22_m);

	// Scale:
	// =====
	double A = params->geometry.A_0;
	IONS->n_m   *= (1/A)*IONS->NCP/params->mesh.DX;
	IONS->nv_m  *= (1/A)*IONS->NCP/params->mesh.DX;
	IONS->P11_m *= (1/A)*IONS->NCP/params->mesh.DX;
	IONS->P22_m *= (1/A)*IONS->NCP/params->mesh.DX;

	// Add finite number for density to avoid zero:
	IONS->n_m += 1E14*CS->volume;
}

void PIC_TYP::calculateDerivedIonMoments(const params_TYP * params, CS_TYP * CS, ionSpecies_TYP * IONS)
{
	double Ma(IONS->M);

	// Ion pressures:
	arma::vec Ppar = IONS->P11_m - (Ma*IONS->nv_m % IONS->nv_m/IONS->n_m);
	arma::vec Pper = IONS->P22_m; // We have neglected perp drift kinetic energy

	// Ion temperatures:
	IONS->Tpar_m = Ppar/(F_E_DS*IONS->n_m);
	IONS->Tper_m = Pper/(F_E_DS*IONS->n_m);
}
