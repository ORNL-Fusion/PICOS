#include <math.h>
#include "particleBC.h"

using namespace std;

particleBC_TYP::particleBC_TYP()
{
    N1_dot = 0;
    N2_dot = 0;
    N5_dot = 0;
    E1_dot = 0;
    E2_dot = 0;
    E5_dot = 0;
    DT     = 0;
}

// =============================================================================
void particleBC_TYP::MPI_AllreduceDouble(const params_TYP * params, double * v)
{
    double recvbuf = 0;

    MPI_Allreduce(v, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, params->mpi.COMM);

    *v = recvbuf;
}

// =============================================================================
void particleBC_TYP::checkBoundaryAndFlag(const params_TYP * params,const CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        for (int aa=0; aa < IONS->size(); aa++)
        {
            // Number of particles is "aa" species:
            // ===================================
            int NSP = IONS->at(aa).NSP;

            // Ion mass:
            // =========
            double Ma = IONS->at(aa).M;

            // Particle loop:
            // ==================================
            #pragma omp parallel for default(none) shared(params, IONS, aa, CS, std::cout) firstprivate(NSP,Ma)
            for(int ii=0; ii<NSP; ii++)
            {
                // left boundary:
                if (IONS->at(aa).X_p(ii) <= params->geometry.LX_min)
                {
                    // Particle flag:
                    IONS->at(aa).f1(ii) = 1;

                    // Particle kinetic energy:
                    double KE = 0.5*Ma*dot(IONS->at(aa).V_p.row(ii), IONS->at(aa).V_p.row(ii));
                    IONS->at(aa).dE1(ii) = KE;
                }

                // Right boundary:
                if (IONS->at(aa).X_p(ii) >= params->geometry.LX_max)
                {
                    // Particle flag:
                    IONS->at(aa).f2(ii) = 1;

                    // Particle kinetic energy:
                    double KE = 0.5*Ma*dot(IONS->at(aa).V_p.row(ii), IONS->at(aa).V_p.row(ii));
                    IONS->at(aa).dE2(ii) = KE;

                }
            } // Particle loop

        } // species loop

    } // Particle MPI guard

}

// =============================================================================
template <typename vec_TYP>
void particleBC_TYP::MPI_OMP_AllreduceVec(const params_TYP * params, vec_TYP * V, double * S)
{
    // Clear S:
    // =======
    (*S) = 0;

    // Number of computational particles per process:
    // ==============================================
    int NSP = V->n_elem;

    // Reduce S over all threads in a single MPI process:
    // ==================================================
    #pragma omp parallel default(none) shared(V, S, std::cout) firstprivate(NSP)
    {
        // Computational particle accumulators:
        double S_private = 0;

        #pragma omp for
        for(int ii=0; ii<NSP; ii++)
        {
            S_private += (*V)(ii);
        }

        #pragma omp critical
        (*S) += S_private;

    } // pragma omp parallel

    // AllReduce S over all PARTICLE MPIs:
    // ==================================
    MPI_AllreduceDouble(params,S);
}

// =============================================================================
void particleBC_TYP::calculateParticleWeight(const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{

    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        // Simulation time step:
        DT = params->DT;

        // Iterate over all ion species:
        // =============================
        for(int ss=0;ss<IONS->size();ss++)
        {
            if (IONS->at(ss).p_BC.BC_type == 1 || IONS->at(ss).p_BC.BC_type == 2)
            {
                // Super particle conversion factor:
                double alpha = IONS->at(ss).NCP;

                // Fueling rate:
                double G = IONS->at(ss).p_BC.G;

                // Computational particle leak rate over all PARTICLE MPIs:
                // ========================================================
                // Store total number of comp particles leaked over all threads and ranks:
                double S1_global = 0;
                double S2_global = 0;

                MPI_OMP_AllreduceVec(params,&IONS->at(ss).f1, &S1_global);
                MPI_OMP_AllreduceVec(params,&IONS->at(ss).f2, &S2_global);

                // Accumulate computational particles and fueling rate:
                // ==================================
                IONS->at(ss).p_BC.S1   += S1_global;
                IONS->at(ss).p_BC.S2   += S2_global;
                IONS->at(ss).p_BC.GSUM += G;

                // Calculate particle weight:
                // ==========================
                // Minimum number of computational particles to trigger fueling:
                double S_min = 3;

                // Total number of computational particles that have leaked:
                double S_total  = IONS->at(ss).p_BC.S1 + IONS->at(ss).p_BC.S2;

                // Check if enough particles have left the domain:
                if ( S_total >= S_min )
                {
                    // Calculate computational particle leak rate:
                    double uN_total = (alpha/DT)*S_total;

                    // Calculate particle weight:
                    double GSUM  = IONS->at(ss).p_BC.GSUM;
                    double a_p_new = GSUM/uN_total;

                    if (a_p_new > 1000)
                    {
                        if (params->mpi.IS_PARTICLES_ROOT)
                        {
                            cout << "S_total:" << S_total << endl;
                            cout << "uN_total:" << uN_total << endl;
                            cout << "a_p_new:" << a_p_new << endl;
                            cout << "GSUM:" << GSUM << endl;
                        }
                        a_p_new = 1000;
                    }
                    IONS->at(ss).p_BC.a_p_new = a_p_new;

                    // Reset accumulators:
                    IONS->at(ss).p_BC.S1   = 0;
                    IONS->at(ss).p_BC.S2   = 0;
                    IONS->at(ss).p_BC.GSUM = 0;
                }

            } // BC

        } //  Species

    } // Particle MPIs

}

void particleBC_TYP::getFluxesAcrossBoundaries(const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        N1_dot = 0;
        N2_dot = 0;
        E1_dot = 0;
        E2_dot = 0;
        DT     = params->DT;

        for (int ss=0; ss<IONS->size();ss++)
        {
            double NCP = IONS->at(ss).NCP;
            int NSP    = IONS->at(ss).NSP;

            #pragma omp parallel default(none) shared(N1_dot, E1_dot, N2_dot, E2_dot, params, IONS, ss, CS, std::cout) firstprivate(NSP,NCP,DT)
            {
                // Private variables:
                double N1_dot_private = 0;
                double N2_dot_private = 0;
                double E1_dot_private = 0;
                double E2_dot_private = 0;

                #pragma omp for
                for(int ii=0; ii<NSP; ii++)
                {
                    // Boundary 1:
                    if ( IONS->at(ss).f1(ii) == 1 )
                    {
                        double dE1 = IONS->at(ss).dE1(ii);
                        double a_p = IONS->at(ss).a_p(ii);

                        // Accumulate fluxes:
                        N1_dot_private += (NCP/DT)*a_p;
                        E1_dot_private += (NCP/DT)*a_p*dE1;

                    }

                    // Boundary 2:
                    if ( IONS->at(ss).f2(ii) == 1 )
                    {
                        double dE2 = IONS->at(ss).dE2(ii);
                        double a_p = IONS->at(ss).a_p(ii);

                        // Accumulate fluxes:
                        N2_dot_private += (NCP/DT)*a_p;
                        E2_dot_private += (NCP/DT)*a_p*dE2;

                    } // if

                } // omp for

                #pragma omp critical
                {
                    N1_dot += N1_dot_private;
                    N2_dot += N2_dot_private;
                    E1_dot += E1_dot_private;
                    E2_dot += E2_dot_private;
                }

            } // omp parallel

        } // Species

        // Reduce over all MPI process
        MPI_AllreduceDouble(params,&N1_dot);
        MPI_AllreduceDouble(params,&N2_dot);
        MPI_AllreduceDouble(params,&E1_dot);
        MPI_AllreduceDouble(params,&E2_dot);

    } // Particle MPI
}

// =============================================================================
void particleBC_TYP::applyParticleReinjection(const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    // Check boundaries:
    // ================
    checkBoundaryAndFlag(params,CS,fields,IONS);

    // Calculate particle and energy flux accross boundaries:
    // ======================================================
    getFluxesAcrossBoundaries(params,CS,fields,IONS);

    // Calculate new particle weight:
    // =============================
    calculateParticleWeight(params,CS,fields,IONS);

    // Apply re-injection:
    // ===================
    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        for(int ss=0;ss<IONS->size();ss++)
        {

                // Number of computational particles per process:
                int NSP(IONS->at(ss).NSP);

                // Ion mass:
                double Ma = IONS->at(ss).M;

    			#pragma omp parallel for default(none) shared(params, CS, fields, IONS, std::cout) firstprivate(NSP, ss, Ma)
                    for(int ii=0; ii<NSP; ii++)
    				{
    					if ( IONS->at(ss).f1(ii) == 1 || IONS->at(ss).f2(ii) == 1 )
    					{

    						// Re-inject particle:
                            // ===================
    						particleReinjection(ii, params, CS, fields,&IONS->at(ss));

                            // Newly injected flag:
                            // ====================
                            IONS->at(ss).f5(ii)  = 1;
                            IONS->at(ss).dE5(ii) = 0.5*Ma*dot(IONS->at(ss).V_p.row(ii), IONS->at(ss).V_p.row(ii));

                            // Reset injection flag:
                            // =====================
    						IONS->at(ss).f1(ii) = 0;
    						IONS->at(ss).f2(ii) = 0;

                            // Reset Exit energy:
                            // =====================
                            IONS->at(ss).dE1(ii) = 0;
                            IONS->at(ss).dE2(ii) = 0;

    					} // flag guard

    				} // pragma omp parallel for

    	} //  Species

	} // Particle MPIs

    // Calculate actual fueling rate and power:
    // =======================================================
    getParticleInjectionRates(params,CS,fields,IONS);

}

void particleBC_TYP::getParticleInjectionRates(const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        N5_dot = 0;
        E5_dot = 0;
        DT     = params->DT;

        for (int ss=0; ss<IONS->size();ss++)
        {
            double NCP = IONS->at(ss).NCP;
            int NSP    = IONS->at(ss).NSP;

            #pragma omp parallel default(none) shared(N5_dot, E5_dot, params, IONS, ss, CS, std::cout) firstprivate(NSP,NCP,DT)
            {
                // Private variables:
                double N5_dot_private = 0;
                double E5_dot_private = 0;

                #pragma omp for
                for(int ii=0; ii<NSP; ii++)
                {
                    // Boundary 1:
                    if ( IONS->at(ss).f5(ii) == 1 )
                    {
                        double dE5 = IONS->at(ss).dE5(ii);
                        double a_p = IONS->at(ss).a_p(ii);

                        // Accumulate fluxes:
                        N5_dot_private += (NCP/DT)*a_p;
                        E5_dot_private += (NCP/DT)*a_p*dE5;

                        // Clear flags:
                        IONS->at(ss).f5(ii)  = 0;
                        IONS->at(ss).dE5(ii) = 0;
                    }

                } // omp for

                #pragma omp critical
                {
                    N5_dot += N5_dot_private;
                    E5_dot += E5_dot_private;
                }

            } // omp parallel

        } // Species

        // Reduce over all MPI process
        MPI_AllreduceDouble(params,&N5_dot);
        MPI_AllreduceDouble(params,&E5_dot);

    } // Particle MPI
}

// =============================================================================
void particleBC_TYP::particleReinjection(int ii, const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, ionSpecies_TYP * IONS)
{
    // Particle velocity:
	// =========================================================================

    // 1: Warm plasma source, 2: NBI or 4: simple-reinjection:
	if (IONS->p_BC.BC_type == 1 || IONS->p_BC.BC_type == 2 || IONS->p_BC.BC_type == 4)
	{
		double T;
		double E;

		if (IONS->p_BC.BC_type == 1 || IONS->p_BC.BC_type == 4) // Warm plasma source
		{
			T = IONS->p_BC.T;
			E = 0;
		}
		if (IONS->p_BC.BC_type == 2) // NBI
		{
			T = IONS->p_BC.T;
			E = IONS->p_BC.E;
		}

		// Mass of ion:
		double Ma = IONS->M;

		// Thermal velocity of source:
		double vT = sqrt(2*F_E_DS*T/Ma);

		// Pitch angle of source:
		double xip = cos(IONS->p_BC.eta);

		// Drift velocity of source:
		double U  = sqrt(2*F_E_DS*E/Ma);
		double Ux = U*xip;
		double Uy = U*sqrt(1 - pow(xip,2));
		double Uz = 0;

		// Thermal spread:
		double sigma_v = vT/sqrt(2);

		// Random number generator:
		std::default_random_engine gen( (params->mpi.MPI_DOMAIN_NUMBER + 1 + time(NULL))*1000 );
		std::uniform_real_distribution<double> Rm(0.0, 1.0);

		// Box muller:
		double R_1 = sigma_v*sqrt(-2*log(Rm(gen)));
		double t_2 = 2*M_PI*Rm(gen);
		double R_3 = sigma_v*sqrt(-2*log(Rm(gen)));
		double t_4 = 2*M_PI*Rm(gen);

		// Thermal component:
		double wx = R_3*cos(t_4);
		double wy = R_1*cos(t_2);
		double wz = R_1*sin(t_2);

		// Total velocity components:
        double v_par = Ux + wx;
        double v_y   = Uy + wy;
        double v_z   = Uz + wz;
        double v_per = sqrt( pow(v_y,2) + pow(v_z,2) );

        // Assign to IONS:
        IONS->V_p(ii,0) = v_par;
        IONS->V_p(ii,1) = v_per;
	}

    // Periodic:
    if (IONS->p_BC.BC_type == 3 )
    {
        // Do nothing
    }

	// Particle position:
    // =========================================================================
    // 1: Warm plasma source, 2: NBI or 4: simple-reinjection:
	if (IONS->p_BC.BC_type == 1 || IONS->p_BC.BC_type == 2 || IONS->p_BC.BC_type == 4)
	{
		// Variables for random number generator:
		arma::vec R = randu(1);
		arma_rng::set_seed_random();
		arma::vec phi = 2.0*M_PI*randu<vec>(1);

		// Gaussian distribution in space:
		double mean_x  = IONS->p_BC.mean_x;
		double sigma_x = IONS->p_BC.sigma_x;
		double new_x   = mean_x;
        new_x += (sigma_x)*sqrt( -2*log(R(0)) )*cos(phi(0));

        // Domain boundaries:
        double LX_min = params->geometry.LX_min;
        double LX_max = params->geometry.LX_max;

        // Correction to prevent injecting out of bounds:
        while( (new_x < LX_min) || (new_x > LX_max) )
		{
			 cout<<"Out of bound, Xp(ii) = "<< new_x << endl;

             // Random number:
			 R = randu(1);
			 phi = 2.0*M_PI*randu<vec>(1);

             // Assign new postion:
			 new_x = mean_x  + (sigma_x)*sqrt( -2*log(R(0)) )*cos(phi(0));

             if ( (new_x > LX_min) || (new_x < LX_max) )
             {
                 cout<< "Corrected X = " << new_x <<  endl;
             }
		}

		IONS->X_p(ii) = new_x;
	}

    //4: Periodic:
	if (IONS->p_BC.BC_type == 3)
	{
        if (IONS->X_p(ii) > params->geometry.LX_max)
        {
            IONS->X_p(ii) = params->geometry.LX_min;
        }

        if (IONS->X_p(ii) < params->geometry.LX_min)
        {
            IONS->X_p(ii) = params->geometry.LX_max;
        }

	}

	// Particle weight:
    // =========================================================================
    // 1: Warm plasma source or 2: NBI
	if (IONS->p_BC.BC_type == 1 || IONS->p_BC.BC_type == 2)
	{
		IONS->a_p(ii) = IONS->p_BC.a_p_new;
	}

    // 3: Periodic or 4: basic re-injection
    if (IONS->p_BC.BC_type == 3 || IONS->p_BC.BC_type == 4)
    {
        // Do nothing
        //IONS->a_p(ii) = 1;
    }

}
