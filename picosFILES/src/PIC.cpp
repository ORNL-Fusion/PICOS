#include <algorithm>

#include "PIC.h"

#if 0
void PIC_TYP::MPI_AllreduceVec(const params_TYP * params, arma::vec * v)
{
	arma::vec recvbuf = zeros(v->n_elem);

	MPI_Allreduce(v->memptr(), recvbuf.memptr(), v->n_elem, MPI_DOUBLE, MPI_SUM, params->mpi.MPI_TOPO);

	*v = recvbuf;
}
#endif

void PIC_TYP::MPI_SendVec(const params_TYP &params, arma::vec &v) const
{
	// Send vector from PARTICLE ROOT to FIELDS ROOT:
	if (params.mpi.IS_PARTICLES_ROOT)
    {
		MPI_Send(v.memptr(), v.n_elem, MPI_DOUBLE, params.mpi.FIELDS_ROOT_WORLD_RANK, PARTICLES_TAG, MPI_COMM_WORLD);
	}

	if (params.mpi.IS_FIELDS_ROOT)
    {
		MPI_Recv(v.memptr(), v.n_elem, MPI_DOUBLE, params.mpi.PARTICLES_ROOT_WORLD_RANK, PARTICLES_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// Broadcast from FIELDS ROOT to all other FIELDS MPIs:
	if (params.mpi.COMM_COLOR == FIELDS_MPI_COLOR)
    {
        MPI_Bcast(v.memptr(), v.n_elem, MPI_DOUBLE, 0, params.mpi.COMM);
    }
}


void PIC_TYP::MPI_ReduceVec(const params_TYP &params, arma::vec &v) const
{
    // Reduce the vector at the PARTICLE ROOT:
    // ======================================
    if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        if (params.mpi.IS_PARTICLES_ROOT)
        {
            MPI_Reduce(MPI_IN_PLACE, v.memptr(), v.n_elem, MPI_DOUBLE, MPI_SUM, 0, params.mpi.COMM);
        }
        else
        {
            // Create receive buffer:
            // ======================
            //arma::vec recvbuf = zeros(v.n_elem);

            MPI_Reduce(v.memptr(), NULL, v.n_elem, MPI_DOUBLE, MPI_SUM, 0, params.mpi.COMM);
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


void PIC_TYP::MPI_Recvvec(const params_TYP &params, arma::vec &field) const
{
	// We send the vector from root process of fields to root process of particles
	arma::vec recvbuf(params.mesh.NX_IN_SIM);
	arma::vec sendbuf(params.mesh.NX_IN_SIM);

	sendbuf = field.subvec(1, params.mesh.NX_IN_SIM);

	if (params.mpi.IS_FIELDS_ROOT)
    {
		MPI_Send(sendbuf.memptr(), params.mesh.NX_IN_SIM, MPI_DOUBLE, params.mpi.PARTICLES_ROOT_WORLD_RANK, 0, MPI_COMM_WORLD);
	}

	if (params.mpi.IS_PARTICLES_ROOT)
    {
		MPI_Recv(recvbuf.memptr(), params.mesh.NX_IN_SIM, MPI_DOUBLE, params.mpi.FIELDS_ROOT_WORLD_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		field.subvec(1, params.mesh.NX_IN_SIM) = recvbuf;
	}

	// Then, the fields is broadcasted to all processes in the particles communicator COMM
	if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
		sendbuf = field.subvec(1, params.mesh.NX_IN_SIM);

		MPI_Bcast(sendbuf.memptr(), params.mesh.NX_IN_SIM, MPI_DOUBLE, 0, params.mpi.COMM);

		field.subvec(1, params.mesh.NX_IN_SIM) = sendbuf;
	}
}

void PIC_TYP::MPI_Recv_AllFields(const params_TYP &params, fields_TYP &fields) const
{
	// Send field data from FIELDS ranks and recieve fields data at PARTICLE ranks
	MPI_Recvvec(params,fields.EX_m);
	MPI_Recvvec(params,fields.BX_m);
	MPI_Recvvec(params,fields.dBX_m);

	if (params.SW.RFheating == 1)
	{
		MPI_Recvvec(params,fields.ddBX_m);
	}
}

void PIC_TYP::fillGhosts(arma::vec &C) const
{
	const int NX = C.n_elem;

	C(0)    = C(1);
	C(NX-1) = C(NX-2);
}

void PIC_TYP::fillGhost_AllFields(const params_TYP &params, fields_TYP &fields) const
{
	fillGhosts(fields.EX_m);
	fillGhosts(fields.BX_m);
	fillGhosts(fields.dBX_m);

	if (params.SW.RFheating == 1)
	{
		fillGhosts(fields.ddBX_m);
	}
}

void PIC_TYP::fill4Ghosts(arma::vec &v) const
{
	const int NX = v.n_elem;

	v.subvec(0,1)       = v.subvec(2,3);
	v.subvec(NX-2,NX-1) = v.subvec(NX-4,NX-3);
}

void PIC_TYP::smooth(arma::vec &v, double as) const
{
	const int NX = v.n_elem;

	arma::vec b = zeros(NX);

	const double wc = 0.75; 	// center weight
	const double ws = 0.125;	// sides weight

	//Step 1: Averaging process
	b.subvec(1, NX-2) = v.subvec(1, NX-2);

	fillGhosts(b);

	b.subvec(1, NX-2) = wc*b.subvec(1, NX-2) + ws*b.subvec(2, NX-1) + ws*b.subvec(0, NX-3);

	//Step 2: Averaged weighted variable estimation.
	v.subvec(1, NX-2) = (1.0 - as)*v.subvec(1, NX-2) + as*b.subvec(1, NX-2);
}

// Constructor:
PIC_TYP::PIC_TYP(const params_TYP &params, CS_TYP &CS, fields_TYP &fields, vector<ionSpecies_TYP> &IONS, electrons_TYP &electrons) :
randoms(picos::random::instances<double, uniform, 0.0, 2*numbers::pi_v<double>> (device()))
{
	// Get latest mesh-defined values from FIELDS ranks:
	// =================================================
	MPI_Recv_AllFields(params, fields);

	// Fill the ghost cells in all fields:
	fillGhost_AllFields(params,fields);

	// Assign cell for all particles:
    // =============================
	assignCell_AllSpecies(params,IONS);

	// Interpolate all fields on all species:
	// ======================================
	interpolateFields_AllSpecies(params,IONS,fields);

	// Interpolate electron temperature on all species:
	// ===============================================
	interpolateElectrons_AllSpecies(params,IONS,electrons);

	// Calculate ion moments and populate mesh-defined ion moments:
	// ============================================================
	// Run 3 dummy cycles to load "n" and "nv" at previous time steps:
	extrapolateMoments_AllSpecies(params,CS,fields,IONS);
    extrapolateMoments_AllSpecies(params,CS,fields,IONS);
    extrapolateMoments_AllSpecies(params,CS,fields,IONS);

//	if (params.mpi.IS_PARTICLES_ROOT)
//	{
		//	cout << IONS->at(0).n_m/CS->volume<< endl;
//	}

    // Set up pre, post and method functions.
    switch (params.advanceParticleMethod)
    {
        case 1:
            pre = [](const double EM, const double Ma, const double vper, double &Z0)->void {};
            post = [](const double EM, const double Ma, double &Z1)->void {};
            method = [](const double qa, const double Ma, const std::array<double, 3> &EM, const std::array<double, 3> &ZN, std::array<double, 3> &F)->void
            {
                // Gather fields:
                const double E    = EM[0];
                const double B    = EM[1];
                const double dB   = EM[2];

                // Gather particle states:
                const double vpar = ZN[1];
                const double vper = ZN[2];

                // Output:
                F[0] = vpar;
                F[2] = 0.5*vper*vpar*dB/B;
                F[1] = -F[2] + (qa/Ma)*E;
            };
            break;
        case 2:
            pre = [](const double EM, const double Ma, const double vper, double &Z0)->void
            {
                Z0 = 0.5*Ma*vper*vper/EM;
            };
            post = [](const double EM, const double Ma, double &Z1)->void
            {
                Z1 = sqrt(2*Z1*EM/Ma);
            };
            method = [](const double qa, const double Ma, const std::array<double, 3> &EM, const std::array<double, 3> &ZN, std::array<double, 3> &F)->void
            {
                // Gather fields:
                const double E    = EM[0];
                const double dB   = EM[2];

                // Gather particle states:
                double vpar = ZN[1];
                double mu = ZN[2];

                // Output:
                F[0] = +vpar;
                F[1] = -(mu/Ma)*dB + (qa/Ma)*E;
                F[2] = 0;
            };
            break;
    }
}

void PIC_TYP::interpolateScalarField(const params_TYP &params, ionSpecies_TYP &IONS, const arma::vec &F_m, arma::vec &F_p) const
{
    const int NX = params.mesh.NX_IN_SIM + 4; //Mesh size along the X axis (considering the gosht cell)

    // Allocate memory and initialize to zero:
    arma::vec F = zeros(NX);

    // Fill in vector F with mesh-defined data:
    F.subvec(1,NX-2) = F_m;

    // Take care of ghost cells:
    fill4Ghosts(F);

    const int iie = IONS.NSP;
    #pragma omp parallel for default(none) shared(params, IONS, F_p, F, iie)
    for(int ii=0; ii<iie; ii++)
    {
        const int ix = IONS.mn(ii) + 2;

        F_p(ii)  = IONS.wxl(ii)*F(ix-1);
        F_p(ii) += IONS.wxc(ii)*F(ix);
        F_p(ii) += IONS.wxr(ii)*F(ix+1);

    }// omp parallel for
}

void PIC_TYP::interpolateFields(const params_TYP &params, ionSpecies_TYP &IONS, const fields_TYP &fields) const
{
	// Need to use SW in order to enable/disable ddBX interpolation

	// Interpolate mesh-defined fields into particles:
	interpolateScalarField(params, IONS, fields.EX_m , IONS.EX_p );
	interpolateScalarField(params, IONS, fields.BX_m , IONS.BX_p );
	interpolateScalarField(params, IONS, fields.dBX_m, IONS.dBX_p);

	if (params.SW.RFheating == 1)
	{
		interpolateScalarField(params, IONS, fields.ddBX_m, IONS.ddBX_p);
	}
}

void PIC_TYP::interpolateFields_AllSpecies(const params_TYP &params, vector<ionSpecies_TYP> &IONS, const fields_TYP &fields) const
{
    if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        // Interpolate all fields on all species:
        // =======================
        for (ionSpecies_TYP &ION : IONS)
        {
			// Interpolate mesh-defined fields into ALL particles locations:
			interpolateFields(params, ION, fields);
		}
	}
}

void PIC_TYP::interpolateElectrons(const params_TYP &params, ionSpecies_TYP &IONS, const electrons_TYP &electrons) const
{
	// Interpolate mesh-defined fields into particles:
	interpolateScalarField(params, IONS, electrons.Te_m , IONS.Te_p);
}

void PIC_TYP::interpolateElectrons_AllSpecies(const params_TYP &params, vector<ionSpecies_TYP> &IONS, const electrons_TYP &electrons) const
{
    if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
	// Interpolate all fields on all species:
	// =======================
        for (ionSpecies_TYP &ION : IONS)
        {
			// Interpolate mesh-defined electron temperature into ALL particles locations:
			interpolateElectrons(params, ION, electrons);
		}
	}
}

void PIC_TYP::interpEM(const params_TYP &params, const fields_TYP &fields, const double xp, std::array<double, 3> &EM) const
{
    // Assign cell:
    double xMin = params.geometry.LX_min;
    double DX   = params.mesh.DX;
    // During RK4, some projections can be go out of bound
    const int m = min(max(static_cast<int> (round( 0.5 + (xp - xMin)/DX )) - 1, 0), params.mesh.NX_IN_SIM - 1);

    // Distance to nearest grid point:
    const double X = params.mesh.nodesX(m) - xp;

    // Assignment function:
    double xnorm = 1.5 + ((X - DX)/DX);
    const double W0 = 0.5*xnorm*xnorm;    // Left:
    xnorm = X/DX;
    const double W1 = 0.75 - xnorm*xnorm; // Center:
    xnorm = 1.5 - ((X + DX)/DX);
    const double W2 = 0.5*xnorm*xnorm;    // Right:

    // Nearest grid point:
    int ix = m + 1;

    // Interpolate:
    // EX:
    double f0 = fields.EX_m(ix - 1);
    double f1 = fields.EX_m(ix);
    double f2 = fields.EX_m(ix + 1);
    EM[0] = f0*W0 + f1*W1 + f2*W2;

    // BX:
    f0 = fields.BX_m(ix - 1);
    f1 = fields.BX_m(ix);
    f2 = fields.BX_m(ix + 1);
    EM[1] = f0*W0 + f1*W1 + f2*W2;

    // dBX:
    f0 = fields.dBX_m(ix - 1);
    f1 = fields.dBX_m(ix);
    f2 = fields.dBX_m(ix + 1);
    EM[2] = f0*W0 + f1*W1 + f2*W2;

}

void PIC_TYP::calculateF(const params_TYP &params, const ionSpecies_TYP &IONS, const std::array<double, 3> &ZN, const std::array<double, 3> &EM, std::array<double, 3> &F) const
{
    // Ion parameters:
    const double qa = IONS.Q;
    const double Ma = IONS.M;

    method(qa, Ma, EM, ZN, F);
}

void PIC_TYP::advanceParticles(const params_TYP &params, fields_TYP &fields, vector<ionSpecies_TYP> &IONS) const
{
    // Get latest mesh-defined values from FIELDS MPIs:
	MPI_Recv_AllFields(params, fields);

	// Fill the ghost cells in all fields:
	fillGhost_AllFields(params, fields);

	// Iterate over all the ion species:
    if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
	{
        // Time step:
        const double DT = params.DT;

        for(auto &ion: IONS)
		{
			// Number of particles:
			const int NSP = ion.NSP;

			// Ion mass:
			const double Ma = ion.M;

			#pragma omp parallel for default(none) shared(params, fields, DT, ion, Ma, std::cout) firstprivate(NSP, F_C_DS)
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
                std::array<double, 3> EM;
                std::array<double, 3> F;
                
                // Extract particle states:
                double x    = ion.X_p(ii);
                double vpar = ion.V_p(ii,0);
                double vper = ion.V_p(ii,1);
                
                // Initialize initial particle state Z0:
                std::array<double, 3> Z0 = {x, vpar, vper};
                
                if ( isnan(Z0[0]) || isnan(Z0[1]) || isnan(Z0[2]) )
                {
                    cout << "*Z0(0), 1 = " << Z0[0] << endl;
                    cout << "*Z0(1), 1 = " << Z0[1] << endl;
                    cout << "*Z0(2), 1 = " << Z0[2] << endl;
                }
                
                // Interpolate fields at current particle postion:
                interpEM(params, fields, Z0[0], EM);
                
                // Select solution method:
                pre(EM[1], Ma, vper, Z0[2]);
                
                // Step 1:
                calculateF(params, ion, Z0, EM, F);
                const std::array<double, 3> dZ1 = {
                    F[0]*DT,
                    F[1]*DT,
                    F[2]*DT
                };
                
                // Step 2:
                std::array<double, 3> ZN = {
                    Z0[0] + dZ1[0]/2,
                    Z0[1] + dZ1[1]/2,
                    Z0[2] + dZ1[2]/2
                };
                interpEM(params, fields, ZN[0], EM);
                calculateF(params, ion, ZN, EM, F);
                const std::array<double, 3> dZ2 = {
                    F[0]*DT,
                    F[1]*DT,
                    F[2]*DT
                };

                // Step 3:
                ZN[0] = Z0[0] + dZ2[0]/2;
                ZN[1] = Z0[1] + dZ2[1]/2;
                ZN[2] = Z0[2] + dZ2[2]/2;;
                interpEM(params, fields, ZN[0], EM);
                calculateF(params, ion, ZN, EM, F);
                const std::array<double, 3> dZ3 = {
                    F[0]*DT,
                    F[1]*DT,
                    F[2]*DT
                };

                // Step 4:
                ZN[0] = Z0[0] + dZ3[0];
                ZN[1] = Z0[1] + dZ3[1];
                ZN[2] = Z0[2] + dZ3[2];
                interpEM(params, fields, ZN[0], EM);
                calculateF(params, ion, ZN, EM, F);
                const std::array<double, 3> dZ4 = {
                    F[0]*DT,
                    F[1]*DT,
                    F[2]*DT
                };

                // Assemble RK4 solution:
                std::array<double, 3> Z1 = {
                    Z0[0] + (dZ1[0] + 2*(dZ2[0] + dZ3[0]) + dZ4[0])/6,
                    Z0[1] + (dZ1[1] + 2*(dZ2[1] + dZ3[1]) + dZ4[1])/6,
                    Z0[2] + (dZ1[2] + 2*(dZ2[2] + dZ3[2]) + dZ4[2])/6
                };

				if ( isnan(ZN[0]) || isnan(ZN[1]) || isnan(ZN[2]) )
				{
					cout << "Z0(0) = " << Z0[0] << endl;
					cout << "Z1(0) = " << Z1[0] << endl;
					cout << "Z1(1) = " << Z1[1] << endl;
					cout << "Z1(2) = " << Z1[2] << endl;
				}

				// Interpolate fields at new particle position:
				interpEM(params, fields, Z1[0], EM);

				// Assign solution to output vector:
                post(EM[1], Ma, Z1[2]);
                // End of RK solution:
                //==============================================================

                // Update new particle states:
                ion.X_p(ii)   = Z1[0];
                ion.V_p(ii,0) = Z1[1]; // vpar
                ion.V_p(ii,1) = Z1[2]; // vper
                ion.mu_p(ii)  = 0.5*Ma*Z1[2]*Z1[2]/EM[1] ; // mu

			} // End of parallel region
		} //structure to iterate over all the ion species.
	}
}

void PIC_TYP::assignCell(const params_TYP &params, ionSpecies_TYP &ION) const
{
	// Clear assignment function:
    ION.wxc.zeros();
    ION.wxl.zeros();
    ION.wxr.zeros();

    const int iie = ION.NSP;
	#pragma omp parallel for default(none) shared(ION, params, std::cout, iie)
    for(int ii=0; ii<iie; ii++)
    {
		// Calculate nearest grid point:
		const double X_p     = ION.X_p(ii);
        const double X_p_min = params.geometry.LX_min;
        const double DX      = params.mesh.DX;
        // Correct "m" near boundaries if out of bound:
        const int m = min(max(static_cast<int> (round( 0.5 + (X_p - X_p_min)/DX )) - 1, 0), params.mesh.NX_IN_SIM - 1);
        // Assign nearest grid point:
        ION.mn(ii) = m;

		// Distance to nearest grid point:
		const double X = params.mesh.nodesX(m) - X_p;

		// Assignment function:
        const double plus = 1.5 + ((X - DX)/DX);
        const double minus = 1.5 - ((X + DX)/DX);
		ION.wxl(ii) = 0.5*plus*plus; // Left:
		ION.wxc(ii) = 0.75 - X*X/(DX*DX);             // Center:
		ION.wxr(ii) = 0.5*minus*minus; // Right:

	} // parallel omp

}

void PIC_TYP::assignCell_AllSpecies(const params_TYP &params, vector<ionSpecies_TYP> &IONS) const
{
    // Assign cell and calculate partial ion moments:
    if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
    // Iterate over all ion species:
    // =============================
        for (ionSpecies_TYP &ION : IONS)
        {
            // Assign cell:
            PIC_TYP::assignCell(params, ION);
        }
	}
}

void PIC_TYP::extrapolateMoments_AllSpecies(const params_TYP &params, CS_TYP &CS, fields_TYP &fields, vector<ionSpecies_TYP> &IONS) const
{
    const int filters = params.filtersPerIterationIons;

	// Iterate over all ion species:
    // =============================
    for (ionSpecies_TYP &ion : IONS)
    {

        // Assign cell and calculate partial ion moments:
        if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
        {
			//Calculate partial moments:
			calculateIonMoments(params,CS,fields,ion);

			// Reduce IONS moments to PARTICLE ROOT:
			// =====================================
			MPI_ReduceVec(params, ion.n_m);
			MPI_ReduceVec(params, ion.nv_m);
			MPI_ReduceVec(params, ion.P11_m);
			MPI_ReduceVec(params, ion.P22_m);

			// Broadcast to all PARTICLE ranks:
			// ================================
			MPI_Bcast(ion.n_m.memptr()  , ion.n_m.size()  , MPI_DOUBLE, 0, params.mpi.COMM);
			MPI_Bcast(ion.nv_m.memptr() , ion.nv_m.size() , MPI_DOUBLE, 0, params.mpi.COMM);
			MPI_Bcast(ion.P11_m.memptr(), ion.P11_m.size(), MPI_DOUBLE, 0, params.mpi.COMM);
			MPI_Bcast(ion.P22_m.memptr(), ion.P22_m.size(), MPI_DOUBLE, 0, params.mpi.COMM);

			// Apply smoothing:
			// ===============
			for (int jj=0; jj<filters; jj++)
			{
			  smooth(ion.n_m  , params.smoothingParameter);
			  smooth(ion.nv_m , params.smoothingParameter);
			  smooth(ion.P11_m, params.smoothingParameter);
			  smooth(ion.P22_m, params.smoothingParameter);
			}

			// Add finite number to density to avoid zero:
			// ============================================
			//IONS->at(ss).n_m += 1E16*CS->volume;

			// Calculate derived ion moments: Tpar_m, Tper_m:
			// ==============================================
			calculateDerivedIonMoments(params, CS, ion);
        }

		// 0th moment at various time levels are sent to fields processes:
        // =============================================================
		// Ion density:
        MPI_SendVec(params, ion.n_m);
        MPI_SendVec(params, ion.n_m_);
        MPI_SendVec(params, ion.n_m__);
        MPI_SendVec(params, ion.n_m___);
	}
}

void PIC_TYP::calculateIonMoments(const params_TYP &params, CS_TYP &CS, fields_TYP &fields, ionSpecies_TYP &ION) const
{
	// Ion density:
	ION.n_m___ = ION.n_m__;
	ION.n_m__  = ION.n_m_;
	ION.n_m_   = ION.n_m;

	// Ion flux:
	ION.nv_m__ = ION.nv_m_;
	ION.nv_m_  = ION.nv_m;

	// Calculate ion moments:
	eim(params,CS,fields,ION);
}

void PIC_TYP::eim(const params_TYP &params, CS_TYP &CS, fields_TYP &fields, ionSpecies_TYP &ION) const
{
	// Ion mass:
	const double Ma = ION.M;

	// Reference magnetic field:
	const double B0 = params.em_IC.BX; // Maybe use the current value at the reference location

	// Clearing content of ion moments:
	// ===============================
	ION.n_m.zeros();
	ION.nv_m.zeros();
	ION.P11_m.zeros();
	ION.P22_m.zeros();

	#pragma omp parallel default(none) shared(params, ION, B0, Ma)
	{
        uniform_random &randuni = randoms[picos::random::thread()];

		// Create private moments:
		// ======================
		arma::vec n   = zeros(params.mesh.NX_IN_SIM + 4);
		arma::vec nv  = zeros(params.mesh.NX_IN_SIM + 4);
		arma::vec P11 = zeros(params.mesh.NX_IN_SIM + 4);
		arma::vec P22 = zeros(params.mesh.NX_IN_SIM + 4);

		// Assemble moments:
		// =================
                const int iie = ION.NSP;
        #pragma omp parallel for
		for(int ii=0; ii<iie; ii++)
		{
			// Nearest grid point:
			const int ix = ION.mn(ii) + 2;

			// Particle velocity:
			const double vpar = ION.V_p(ii,0);
			const double vper = ION.V_p(ii,1);

			// vx component:
			const double vy = vper*cos(randuni());

			// Particle-defined magnetic field:
			const double B = ION.BX_p(ii);

			// Compression factor:
			//double c = B/B0;
			const double c = 1;

			// Particle weight:
			/*
			if (params->currentTime == 0)
			{
				IONS->a_p(ii) = 1/c;
			}
			*/

			const double ac = ION.a_p(ii);
            const double wl = ION.wxl(ii);
            const double wc = ION.wxc(ii);
            const double wr = ION.wxr(ii);

			// Density:
			n(ix-1) += wl*ac;
			n(ix)   += wc*ac;
			n(ix+1) += wr*ac;

			// Particle flux density:
            const double acvpar = ac*vpar;
			nv(ix-1) += wl*acvpar;
			nv(ix) 	 += wc*acvpar;
			nv(ix+1) += wr*acvpar;

			// Stress tensor P11:
            const double acmvpar2 = acvpar*Ma*vpar;
			P11(ix-1) += wl*acmvpar2;
			P11(ix)   += wc*acmvpar2;
			P11(ix+1) += wr*acmvpar2;

			// Stress tensor P22:
            const double acmvy2 = ac*Ma*vy*vy;
			P22(ix-1) += wl*acmvy2;
			P22(ix)   += wc*acmvy2;
			P22(ix+1) += wr*acmvy2;
		}

		// Reduce partial moments from each thread:
		// ========================================
		#pragma omp critical (update_ion_moments)
		{
			ION.n_m.subvec(1,params.mesh.NX_IN_SIM)   += n.subvec(2,params.mesh.NX_IN_SIM + 1);
			ION.nv_m.subvec(1,params.mesh.NX_IN_SIM)  += nv.subvec(2,params.mesh.NX_IN_SIM + 1);
			ION.P11_m.subvec(1,params.mesh.NX_IN_SIM) += P11.subvec(2,params.mesh.NX_IN_SIM + 1);
			ION.P22_m.subvec(1,params.mesh.NX_IN_SIM) += P22.subvec(2,params.mesh.NX_IN_SIM + 1);
		}

	}//End of the parallel region

	// Ghost contributions:
	// ====================
	fill4Ghosts(ION.n_m);
	fill4Ghosts(ION.nv_m);
	fill4Ghosts(ION.P11_m);
	fill4Ghosts(ION.P22_m);

	// Apply compression factor:
	// ========================
	const arma::vec c = fields.BX_m/B0;
	ION.n_m   = ION.n_m%c;
	ION.nv_m  = ION.nv_m%c;
	ION.P11_m = ION.P11_m%c;
	ION.P22_m = ION.P22_m%c;

	// Scale:
	// =====
	const double A = params.geometry.A_0;
	ION.n_m   *= (1/A)*ION.NCP/params.mesh.DX;
	ION.nv_m  *= (1/A)*ION.NCP/params.mesh.DX;
	ION.P11_m *= (1/A)*ION.NCP/params.mesh.DX;
	ION.P22_m *= (1/A)*ION.NCP/params.mesh.DX;

}

void PIC_TYP::calculateDerivedIonMoments(const params_TYP &params, CS_TYP &CS, ionSpecies_TYP &ION) const
{
	// Ion pressures:

	// Ion temperatures:
	ION.Tpar_m = (ION.P11_m - (ION.M*ION.nv_m % ION.nv_m/ION.n_m))/(F_E_DS*ION.n_m);
	ION.Tper_m = ION.P22_m/(F_E_DS*ION.n_m);
}
