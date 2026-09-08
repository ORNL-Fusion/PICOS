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

    // Spatial increment:
	dx = params->mesh.DX;

    // Electric field:
    EX_m.zeros(NX_S);

    // Electrostatic potential and charge density:
    Phi_m.zeros(NX_T);
    chargeDensity.zeros(NX_T);

    // Reformulated Poisson work arrays:
    ionDensity.zeros(NX_T);
    electronDensity.zeros(NX_T);
    stressDifference.zeros(NX_T);
    divStressDifference.zeros(NX_T);
    plasmaFrequencySquared.zeros(NX_T);
}

// Fill ghost cells:
// ============================================================================================
void fields_solver_TYP::fillGhosts(arma::vec * C)
{
	int NX = C->n_elem;

	(*C)(0)    = (*C)(1);
	(*C)(NX-1) = (*C)(NX-2);
}

void fields_solver_TYP::fillPeriodicGhosts(arma::vec * C)
{
	int NX = C->n_elem;

	(*C)(0)    = (*C)(NX-2);
	(*C)(NX-1) = (*C)(1);
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

void fields_solver_TYP::smoothPeriodic(arma::vec * v, double as)
{
	int NX = v->n_elem;

	arma::vec b = zeros(NX);

	double wc(0.5);
	double ws(0.25);

	b.subvec(1, NX-2) = v->subvec(1, NX-2);
	fillPeriodicGhosts(&b);
	b.subvec(1, NX-2) = wc*b.subvec(1, NX-2) + ws*b.subvec(2, NX-1) + ws*b.subvec(0, NX-3);
	v->subvec(1, NX-2) = (1.0 - as)*v->subvec(1, NX-2) + as*b.subvec(1, NX-2);
	fillPeriodicGhosts(v);
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
	if (params->SW.fieldSolveModel == FIELD_SOLVE_REFORMULATED_POISSON)
	{
		advanceEfieldReformulatedPoisson(params, fields, CS, IONS);
	}
	else if (params->SW.fieldSolveModel == FIELD_SOLVE_POISSON)
	{
		advanceEfieldPoisson(params, fields, CS, IONS);
	}
	else
	{
		advanceEfieldOhmLaw(params, fields, CS, IONS, electrons);
	}
}

void fields_solver_TYP::advanceEfieldOhmLaw(const params_TYP * params, fields_TYP * fields, CS_TYP * CS, vector<ionSpecies_TYP> * IONS, electrons_TYP * electrons)
{
	if (params->mpi.COMM_COLOR == FIELDS_MPI_COLOR)
	{
		// Indices of subdomain:
		unsigned int iIndex = params->mpi.iIndex;
		unsigned int fIndex = params->mpi.fIndex;

		// Initialize the electron density prior to the accumulation:
		ne.zeros();
		ne_.zeros();
		ne__.zeros();
		ne___.zeros();

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

		// Electric field based on Ohm's law:
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

void fields_solver_TYP::solveDirichletPoisson(const params_TYP * params, const arma::vec * rho, arma::vec * phi) const
{
	const int N = params->mesh.NX_IN_SIM;

	phi->zeros(N + 2);

	if (N <= 0)
	{
		return;
	}

	const double phiLeft = poissonBoundaryPotential(params, false);
	const double phiRight = poissonBoundaryPotential(params, true);

	(*phi)(0) = phiLeft;
	(*phi)(N + 1) = phiRight;

	arma::vec lower = ones(N);
	arma::vec diag = -2.0*ones(N);
	arma::vec upper = ones(N);
	arma::vec rhs = zeros(N);

	lower(0) = 0.0;
	upper(N - 1) = 0.0;

	for (int ii=0; ii<N; ii++)
	{
		rhs(ii) = -(*rho)(ii + 1)*dx*dx/F_EPSILON_DS;
	}

	rhs(0) -= phiLeft;
	rhs(N - 1) -= phiRight;

	for (int ii=1; ii<N; ii++)
	{
		const double m = lower(ii)/diag(ii - 1);
		diag(ii) -= m*upper(ii - 1);
		rhs(ii) -= m*rhs(ii - 1);
	}

	arma::vec solution = zeros(N);
	solution(N - 1) = rhs(N - 1)/diag(N - 1);

	for (int ii=N - 2; ii>=0; ii--)
	{
		solution(ii) = (rhs(ii) - upper(ii)*solution(ii + 1))/diag(ii);
	}

	phi->subvec(1,N) = solution;
}

void fields_solver_TYP::solvePeriodicPoisson(const params_TYP * params, const arma::vec * rho, arma::vec * phi) const
{
	const int N = params->mesh.NX_IN_SIM;

	phi->zeros(N + 2);

	if (N <= 1)
	{
		return;
	}

	arma::vec rhs = -rho->subvec(1,N)*dx*dx/F_EPSILON_DS;
	rhs -= mean(rhs);

	arma::cx_vec rhsHat = fft(arma::cx_vec(rhs, zeros(N)));
	arma::cx_vec phiHat = zeros<cx_vec>(N);

	for (int kk=1; kk<N; kk++)
	{
		const double lambda = -4.0*pow(sin(M_PI*(double)kk/(double)N), 2.0);
		phiHat(kk) = rhsHat(kk)/lambda;
	}

	arma::cx_vec solution = ifft(phiHat);
	phi->subvec(1,N) = real(solution);

	(*phi)(0) = (*phi)(N);
	(*phi)(N + 1) = (*phi)(1);
}

double fields_solver_TYP::poissonBoundaryPotential(const params_TYP * params, bool rightBoundary) const
{
	double phi = rightBoundary ? params->em_IC.phiRight : params->em_IC.phiLeft;

	if (params->em_IC.poissonBCModel == POISSON_BC_SHEATH)
	{
		phi -= params->em_IC.sheathCoefficient*params->f_IC.Te/F_E_DS;
	}

	return phi;
}

void fields_solver_TYP::advanceEfieldPoisson(const params_TYP * params, fields_TYP * fields, CS_TYP * CS, vector<ionSpecies_TYP> * IONS)
{
	if (params->mpi.COMM_COLOR == FIELDS_MPI_COLOR)
	{
		// Indices of subdomain:
		unsigned int iIndex = params->mpi.iIndex;
		unsigned int fIndex = params->mpi.fIndex;

		chargeDensity.zeros();
		fields->EX_m.zeros();
		fields->Phi_m.zeros();

		// rho is normalized by n0*q0. The mesh moments store physical density
		// times CS->volume, so divide by CS->density*CS->volume.
		for(int ss=0; ss<params->numberOfParticleSpecies; ss++)
		{
			chargeDensity.subvec(iIndex, fIndex) += IONS->at(ss).Q*IONS->at(ss).n_m.subvec(iIndex, fIndex)/(CS->density*CS->volume);
		}

		MPI_Allgathervec(params, &chargeDensity);
		if (params->em_IC.poissonBCModel == POISSON_BC_PERIODIC)
		{
			fillPeriodicGhosts(&chargeDensity);
		}
		else
		{
			fillGhosts(&chargeDensity);
		}

		for (int jj=0; jj<params->filtersPerIterationFields; jj++)
		{
			if (params->em_IC.poissonBCModel == POISSON_BC_PERIODIC)
			{
				smoothPeriodic(&chargeDensity, params->smoothingParameter);
			}
			else
			{
				smooth(&chargeDensity, params->smoothingParameter);
			}
		}

		if (params->em_IC.poissonBCModel == POISSON_BC_PERIODIC)
		{
			solvePeriodicPoisson(params, &chargeDensity, &fields->Phi_m);
		}
		else
		{
			solveDirichletPoisson(params, &chargeDensity, &fields->Phi_m);
		}

		for (int ii=1; ii<=params->mesh.NX_IN_SIM; ii++)
		{
			fields->EX_m(ii) = -(fields->Phi_m(ii + 1) - fields->Phi_m(ii - 1))/(2.0*dx);
		}
		if (params->em_IC.poissonBCModel == POISSON_BC_PERIODIC)
		{
			fields->EX_m(0) = fields->EX_m(params->mesh.NX_IN_SIM);
			fields->EX_m(params->mesh.NX_IN_SIM + 1) = fields->EX_m(1);
		}
		else
		{
			fillGhosts(&fields->EX_m);
		}

		for (int jj=0; jj<params->filtersPerIterationFields; jj++)
		{
			if (params->em_IC.poissonBCModel == POISSON_BC_PERIODIC)
			{
				smoothPeriodic(&fields->EX_m, params->smoothingParameter);
			}
			else
			{
				smooth(&fields->EX_m, params->smoothingParameter);
			}
		}

		#ifdef CHECKS_ON
		if(!fields->EX_m.is_finite() || !fields->Phi_m.is_finite())
		{
			cout << "Non finite values in electrostatic Poisson field solve" << endl;
			MPI_Abort(params->mpi.MPI_TOPO, -110);
		}
		#endif
	}

	// Send E and Phi to PARTICLE ranks for pusher and sheath diagnostics.
	MPI_SendVec(params,&fields->EX_m);
	MPI_SendVec(params,&fields->Phi_m);
}

void fields_solver_TYP::advanceEfieldReformulatedPoisson(const params_TYP * params, fields_TYP * fields, CS_TYP * CS, vector<ionSpecies_TYP> * IONS)
{
	if (params->mpi.COMM_COLOR == FIELDS_MPI_COLOR)
	{
		// Indices of subdomain:
		unsigned int iIndex = params->mpi.iIndex;
		unsigned int fIndex = params->mpi.fIndex;

		const arma::vec previousEX = fields->EX_m;

		ionDensity.zeros();
		electronDensity.zeros();
		stressDifference.zeros();
		divStressDifference.zeros();
		plasmaFrequencySquared.zeros();

		double referenceIonMass = 0.0;
		double referenceElectronMass = 0.0;

		for(int ss=0; ss<params->numberOfParticleSpecies; ss++)
		{
			const ionSpecies_TYP &ion = IONS->at(ss);
			const double absZ = fabs(ion.Z);
			if (absZ <= double_zero || ion.M <= double_zero)
			{
				continue;
			}

			if (ion.Z > 0.0 && referenceIonMass <= double_zero)
			{
				referenceIonMass = ion.M;
			}
			if (ion.Z < 0.0 && referenceElectronMass <= double_zero)
			{
				referenceElectronMass = ion.M;
			}

			const arma::vec density = absZ*ion.n_m.subvec(iIndex, fIndex)/(CS->density*CS->volume);
			const arma::vec parallelSecondMoment = absZ*ion.P11_m.subvec(iIndex, fIndex)/(ion.M*CS->density*CS->volume);

			if (ion.Z > 0.0)
			{
				ionDensity.subvec(iIndex, fIndex) += density;
				stressDifference.subvec(iIndex, fIndex) += parallelSecondMoment;
			}
			else if (ion.Z < 0.0)
			{
				electronDensity.subvec(iIndex, fIndex) += density;
				stressDifference.subvec(iIndex, fIndex) -= parallelSecondMoment;
			}
		}

		MPI_Allgathervec(params, &ionDensity);
		MPI_Allgathervec(params, &electronDensity);
		MPI_Allgathervec(params, &stressDifference);

		if (params->em_IC.poissonBCModel == POISSON_BC_PERIODIC)
		{
			fillPeriodicGhosts(&ionDensity);
			fillPeriodicGhosts(&electronDensity);
			fillPeriodicGhosts(&stressDifference);
		}
		else
		{
			fillGhosts(&ionDensity);
			fillGhosts(&electronDensity);
			fillGhosts(&stressDifference);
		}

		for (int jj=0; jj<params->filtersPerIterationFields; jj++)
		{
			if (params->em_IC.poissonBCModel == POISSON_BC_PERIODIC)
			{
				smoothPeriodic(&ionDensity, params->smoothingParameter);
				smoothPeriodic(&electronDensity, params->smoothingParameter);
				smoothPeriodic(&stressDifference, params->smoothingParameter);
			}
			else
			{
				smooth(&ionDensity, params->smoothingParameter);
				smooth(&electronDensity, params->smoothingParameter);
				smooth(&stressDifference, params->smoothingParameter);
			}
		}

		for (int ii=1; ii<=params->mesh.NX_IN_SIM; ii++)
		{
			divStressDifference(ii) = (stressDifference(ii + 1) - stressDifference(ii - 1))/(2.0*dx);
		}

		if (params->em_IC.poissonBCModel == POISSON_BC_PERIODIC)
		{
			divStressDifference(0) = divStressDifference(params->mesh.NX_IN_SIM);
			divStressDifference(params->mesh.NX_IN_SIM + 1) = divStressDifference(1);
		}
		else
		{
			fillGhosts(&divStressDifference);
		}

		if (referenceIonMass <= double_zero || referenceElectronMass <= double_zero)
		{
			cout << "Reformulated Poisson requires at least one positive-Z and one negative-Z self-consistent species." << endl;
			MPI_Abort(params->mpi.MPI_TOPO, -110);
		}

		const double epsilon = max(referenceElectronMass/referenceIonMass, double_zero);
		const double lambda = (params->em_IC.reformulatedPoissonLambda > 0.0) ?
			params->em_IC.reformulatedPoissonLambda : sqrt(max(F_EPSILON_DS, double_zero));
		const double lambda2 = lambda*lambda;
		const bool quasiNeutral = (params->em_IC.reformulatedPoissonQuasiNeutral != 0) || (lambda2 <= double_zero);

		fields->EX_m.zeros();
		fields->Phi_m.zeros();

		for (int ii=1; ii<=params->mesh.NX_IN_SIM; ii++)
		{
			const double densityDenominator = max(ionDensity(ii) + electronDensity(ii)/epsilon, double_zero);

			if (quasiNeutral)
			{
				fields->EX_m(ii) = divStressDifference(ii)/densityDenominator;
			}
			else
			{
				plasmaFrequencySquared(ii) = densityDenominator/lambda2;
				const double rhs = divStressDifference(ii)/lambda2;
				fields->EX_m(ii) = (previousEX(ii) + params->DT*rhs)/(1.0 + params->DT*plasmaFrequencySquared(ii));
			}
		}

		if (params->em_IC.poissonBCModel == POISSON_BC_PERIODIC)
		{
			fields->EX_m(0) = fields->EX_m(params->mesh.NX_IN_SIM);
			fields->EX_m(params->mesh.NX_IN_SIM + 1) = fields->EX_m(1);
		}
		else
		{
			fillGhosts(&fields->EX_m);
		}

		for (int jj=0; jj<params->filtersPerIterationFields; jj++)
		{
			if (params->em_IC.poissonBCModel == POISSON_BC_PERIODIC)
			{
				smoothPeriodic(&fields->EX_m, params->smoothingParameter);
			}
			else
			{
				smooth(&fields->EX_m, params->smoothingParameter);
			}
		}

		if (params->em_IC.poissonBCModel == POISSON_BC_PERIODIC)
		{
			fields->EX_m.subvec(1,params->mesh.NX_IN_SIM) -= mean(fields->EX_m.subvec(1,params->mesh.NX_IN_SIM));
			fields->EX_m(0) = fields->EX_m(params->mesh.NX_IN_SIM);
			fields->EX_m(params->mesh.NX_IN_SIM + 1) = fields->EX_m(1);
		}

		if (params->em_IC.poissonBCModel == POISSON_BC_PERIODIC)
		{
			fields->Phi_m(1) = 0.0;
			for (int ii=2; ii<=params->mesh.NX_IN_SIM; ii++)
			{
				fields->Phi_m(ii) = fields->Phi_m(ii - 1) - 0.5*(fields->EX_m(ii) + fields->EX_m(ii - 1))*dx;
			}
			fields->Phi_m.subvec(1,params->mesh.NX_IN_SIM) -= mean(fields->Phi_m.subvec(1,params->mesh.NX_IN_SIM));
			fields->Phi_m(0) = fields->Phi_m(params->mesh.NX_IN_SIM);
			fields->Phi_m(params->mesh.NX_IN_SIM + 1) = fields->Phi_m(1);
		}
		else
		{
			fields->Phi_m(0) = poissonBoundaryPotential(params, false);
			for (int ii=1; ii<=params->mesh.NX_IN_SIM + 1; ii++)
			{
				fields->Phi_m(ii) = fields->Phi_m(ii - 1) - 0.5*(fields->EX_m(ii) + fields->EX_m(ii - 1))*dx;
			}
		}

		#ifdef CHECKS_ON
		if(!fields->EX_m.is_finite() || !fields->Phi_m.is_finite())
		{
			cout << "Non finite values in reformulated electrostatic Poisson field solve" << endl;
			MPI_Abort(params->mpi.MPI_TOPO, -110);
		}
		#endif
	}

	// Send E and reconstructed Phi to PARTICLE ranks for pusher and sheath diagnostics.
	MPI_SendVec(params,&fields->EX_m);
	MPI_SendVec(params,&fields->Phi_m);
}
