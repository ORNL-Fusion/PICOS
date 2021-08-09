#include "PIC.h"

PIC_TYP::PIC_TYP()
{
}

/*
void PIC_TYP::MPI_AllreduceVec(const params_TYP * params, arma::vec * v)
{
	arma::vec recvbuf = zeros(v->n_elem);

	MPI_Allreduce(v->memptr(), recvbuf.memptr(), v->n_elem, MPI_DOUBLE, MPI_SUM, params->mpi.MPI_TOPO);

	*v = recvbuf;
}
*/

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

void PIC_TYP::fillGhosts(arma::vec * C)
{
	int NX = C->n_elem;

	(*C)(0)    = (*C)(1);
	(*C)(NX-1) = (*C)(NX-2);
}

void PIC_TYP::fill4Ghosts(arma::vec * v)
{
	//int N = v->n_elem;
	//v->subvec(N-2,N-1) = v->subvec(2,3);
	//v->subvec(0,1) = v->subvec(N-4,N-3);
}

void PIC_TYP::interpolateScalarField(const params_TYP * params, ionSpecies_TYP * IONS, arma::vec * F_m, arma::vec * F_p)
{
    // Triangular Shape Cloud (TSC) scheme. See Sec. 5-3-2 of R. Hockney and J. Eastwood, Computer Simulation Using Particles.
    //		wxl		   wxc		wxr
    // --------*------------*--------X---*--------
    //				    0       x

    //wxc = 0.75 - (x/H)^2
    //wxr = 0.5*(1.5 - abs(x)/H)^2
    //wxl = 0.5*(1.5 - abs(x)/H)^2

    int NX =  params->mesh.NX_IN_SIM + 4; //Mesh size along the X axis (considering the gosht cell)
    int NSP(IONS->NSP);

    // Allocate memory and initialize to zero:
    arma::vec F = zeros(NX);

    // Fill in vector F with mesh-defined data:
    F.subvec(1,NX-2) = *F_m;

    // Take care of ghost cells:
    fill4Ghosts(&F);
    // fillGhosts(&F);

    #pragma omp parallel for default(none) shared(params, IONS, F_p, F) firstprivate(NSP)
    for(int ii=0; ii<NSP; ii++)
    {
        int ix = IONS->mn(ii) + 2;

        (*F_p)(ii) += IONS->wxl(ii)*F(ix-1);
        (*F_p)(ii) += IONS->wxc(ii)*F(ix);
        (*F_p)(ii) += IONS->wxr(ii)*F(ix+1);

    }//End of the parallel region
}

void PIC_TYP::interpEM(const params_TYP * params, CS_TYP * CS, const ionSpecies_TYP * IONS, const fields_TYP * fields, arma::rowvec * ZN, arma::rowvec * EM)
{
    // Assign cell:
    double xp   = (*ZN)(0); // arma::as_scalar
    double xMin = 0;
    double dx   = params->mesh.DX;
    int m       = round( 0.5 + (xp - xMin)/dx ) - 1;

    // Distance to nearest grid point:
    double X    = params->mesh.nodesX(m) - xp;

    // Assignment function:
    arma::vec W = zeros<vec>(3);
    // Left:
    W(0) = 0.5*pow( 1.5 + ((X - dx)/dx) ,2);
    // Center:
    W(1) = 0.75 - pow(X/dx,2);
    // Right:
    W(2) = 0.5*pow(1.5 - ((X + dx)/dx) ,2 );

    // Nearest grid point:
    int ix = m + 1;

    // Temporary storage for interpolated fields:
    arma::vec f = zeros<vec>(3);

    // Interpolate:
    // EX:
    f(0) = fields->EX_m(ix - 1);
    f(1) = fields->EX_m(ix - 0);
    f(2) = fields->EX_m(ix + 1);
    (*EM)(0) = arma::dot(f,W);

    // BX:
    f(0) = fields->BX_m(ix - 1);
    f(1) = fields->BX_m(ix - 0);
    f(2) = fields->BX_m(ix + 1);
    (*EM)(1) = arma::dot(f,W);


    int m_test = 10000;
    if (params->mpi.IS_PARTICLES_ROOT)
    {
        if (m == m_test)
        {
            cout<< " " << endl;
            cout << "xp ="<< xp*CS->length << endl;
            cout << "m = "<< m << endl;

            cout << "sum(W) = " << sum(W) << endl;

            cout << "BX = "<< (*EM)(1)*CS->bField << endl;
            cout<< " " << endl;

            cout << "f(0) = "<< f(0)*CS->bField << endl;
            cout << "f(1) = "<< f(1)*CS->bField << endl;
            cout << "f(2) = "<< f(2)*CS->bField << endl;
        }
    }

    // dBX:
    f(0) = fields->dBX_m(ix - 1);
    f(1) = fields->dBX_m(ix - 0);
    f(2) = fields->dBX_m(ix + 1);
    (*EM)(2) = arma::dot(f,W);

}

void PIC_TYP::calculateF(const params_TYP * params, CS_TYP * CS, const ionSpecies_TYP * IONS, arma::rowvec * ZN, arma::rowvec * EM, arma::rowvec * F)
{
    // Ion parameters:
    double qa = IONS->Q;
    double Ma = IONS->M;

    // Gather particle states:
    double vpar = (*ZN)(1);
    double vper = (*ZN)(2);
    double E    = (*EM)(0);
    double B    = (*EM)(1);
    double dB   = (*EM)(2);

    // Output:
    (*F)(0) = vpar;
    (*F)(1) = -0.5*vper*vper*dB/B + qa*E/Ma;
    (*F)(2) = +0.5*vper*vpar*dB/B;
}

void PIC_TYP::advanceParticles(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    // Get latest mesh-defined values from FIELDS MPIs:
    MPI_Recvvec(params,&fields->EX_m);
    MPI_Recvvec(params,&fields->BX_m);
    MPI_Recvvec(params,&fields->dBX_m);

    // Fill the ghost cells:
	fillGhosts(&fields->EX_m);
    fillGhosts(&fields->BX_m);
    fillGhosts(&fields->dBX_m);

	// Iterate over all the ion species:
	for(int ss=0; ss<IONS->size(); ss++)
	{
		if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
		{
            // Allocate memory for field variables:
            arma::vec EX_p  = zeros(IONS->at(ss).NSP, 1);
			arma::vec BX_p  = zeros(IONS->at(ss).NSP, 1);
            arma::vec dBX_p = zeros(IONS->at(ss).NSP, 1);

            // Interpolate mesh-defined fields into particles:
            interpolateScalarField(params, &IONS->at(ss), &fields->EX_m , &EX_p );
            interpolateScalarField(params, &IONS->at(ss), &fields->BX_m , &BX_p );
            interpolateScalarField(params, &IONS->at(ss), &fields->dBX_m, &dBX_p);

			IONS->at(ss).EX_p  = EX_p;
			IONS->at(ss).BX_p  = BX_p;
            IONS->at(ss).dBX_p = dBX_p;

			//Once the electric and magnetic fields have been interpolated to the ions' positions we advance the ions' velocities.
			int NSP = IONS->at(ss).NSP;

            // Time step:
            double DT = params->DT;

            // ******* what is this used for?
            double A = IONS->at(ss).Q*params->DT/IONS->at(ss).M; // A = \alpha in the dimensionless equation for the ions' velocity. (Q*NCP/M*NCP=Q/M)

			#pragma omp parallel for default(none) shared(IONS, params, fields, DT, ss) firstprivate(A, NSP, F_C_DS)
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

                // Extract particle states:
                double x    = IONS->at(ss).X_p(ii);
                double vpar = IONS->at(ss).V_p(ii,0);
                double vper = IONS->at(ss).V_p(ii,1);

                // Extract particle-defined fields:
                double E  = IONS->at(ss).EX_p(ii);
                double B  = IONS->at(ss).BX_p(ii);
                double dB = IONS->at(ss).dBX_p(ii);

                // Assemble vectors to use in RK4 method:
                arma::rowvec Z0  = {x, vpar, vper};
                arma::rowvec ZN  = zeros<rowvec>(3);
                arma::rowvec EM  = {E, B   , dB  };
                arma::rowvec F   = zeros<rowvec>(3);
                arma::rowvec dZ1 = zeros<rowvec>(3);
                arma::rowvec dZ2 = zeros<rowvec>(3);
                arma::rowvec dZ3 = zeros<rowvec>(3);
                arma::rowvec dZ4 = zeros<rowvec>(3);

                // Step 1:
                ZN = Z0;
                interpEM(params, CS, &IONS->at(ss), fields, &ZN, &EM);
                calculateF(params, CS, &IONS->at(ss), &ZN, &EM, &F);
                dZ1 = F*DT;

                // Step 2:
                ZN = Z0 + dZ1/2;
                interpEM(params, CS, &IONS->at(ss), fields, &ZN, &EM);
                calculateF(params, CS, &IONS->at(ss), &ZN, &EM, &F);
                dZ2 = F*DT;

                // Step 3:
                ZN = Z0 + dZ2/2;
                interpEM(params, CS, &IONS->at(ss), fields, &ZN, &EM);
                calculateF(params, CS, &IONS->at(ss), &ZN, &EM, &F);
                dZ3 = F*DT;

                // Step 4:
                ZN = Z0 + dZ3;
                interpEM(params, CS, &IONS->at(ss), fields, &ZN, &EM);
                calculateF(params, CS, &IONS->at(ss), &ZN, &EM, &F);
                dZ4 = F*DT;

                // Assemble RK4 solution:
                ZN = Z0 + (dZ1 + 2*dZ2 + 2*dZ3 + dZ4)/6;

                // End of RK solution:
                //==============================================================

                // Output data:
                IONS->at(ss).X_p(ii)   = ZN(0);
                IONS->at(ss).V_p(ii,0) = ZN(1);
                IONS->at(ss).V_p(ii,1) = ZN(2);

			} // End of parallel region
		}
	}//structure to iterate over all the ion species.
}
