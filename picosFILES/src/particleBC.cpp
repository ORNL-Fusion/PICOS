#include <algorithm>

#include "particleBC.h"

using namespace std;

particleBC_TYP::particleBC_TYP() :
dot_({0, 0, 0, 0, 0, 0}),
randoms_2pi(picos::random::instances<double, uniform, 0.0, 2*numbers::pi_v<double>> (device())),
randoms_one(picos::random::instances<double, uniform, 0.0, 1.0> (device())){}

// =============================================================================
void particleBC_TYP::checkBoundaryAndFlag(const params_TYP &params,const CS_TYP &CS, fields_TYP &fields, vector<ionSpecies_TYP> &IONS) const
{
    if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        for (ionSpecies_TYP &ion : IONS)
        {
            // Ion mass:
            // =========
            const double Ma = ion.M;

            // Particle loop:
            // ==================================
            const int iie = ion.NSP;
            #pragma omp parallel for default(none) shared(params, ion, Ma, iie)
            for(int ii=0; ii<iie; ii++)
            {
                // left boundary:
                if (ion.X_p(ii) <= params.geometry.LX_min)
                {
                    // Particle flag:
                    ion.f1(ii) = 1;

                    // Particle kinetic energy:
                    const double KE = 0.5*Ma*dot(ion.V_p.row(ii), ion.V_p.row(ii));
                    ion.dE1(ii) = KE;
                }

                // Right boundary:
                if (ion.X_p(ii) >= params.geometry.LX_max)
                {
                    // Particle flag:
                    ion.f2(ii) = 1;

                    // Particle kinetic energy:
                    const double KE = 0.5*Ma*dot(ion.V_p.row(ii), ion.V_p.row(ii));
                    ion.dE2(ii) = KE;

                }
            } // Particle loop

        } // species loop

    } // Particle MPI guard

}

// =============================================================================
void particleBC_TYP::calculateParticleWeight(const params_TYP &params, const CS_TYP &CS, fields_TYP &fields, vector<ionSpecies_TYP> &IONS) const
{

    if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        // Simulation time step:
        const double DT = params.DT;

        // Iterate over all ion species:
        // =============================
        for (ionSpecies_TYP &ion : IONS)
        {
            if (ion.p_BC.BC_type == 1 || ion.p_BC.BC_type == 2)
            {
                // Computational particle leak rate over all PARTICLE MPIs:
                // ========================================================
                // Store total number of comp particles leaked over all threads and ranks:
                std::array<double, 2> S_global;

                MPI_OMP_AllreduceVec(params, ion.f1, ion.f2, S_global);

                // Accumulate computational particles and fueling rate:
                // ==================================
                ion.p_BC.S1   += S_global[0];
                ion.p_BC.S2   += S_global[1];
                ion.p_BC.GSUM += ion.p_BC.G;

                // Calculate particle weight:
                // ==========================
                // Minimum number of computational particles to trigger fueling:
                const double S_min = 3;

                // Total number of computational particles that have leaked:
                const double S_total  = ion.p_BC.S1 + ion.p_BC.S2;

                // Check if enough particles have left the domain:
                if ( S_total >= S_min )
                {
                    // Calculate computational particle leak rate:
                    const double uN_total = (ion.NCP/DT)*S_total;

                    // Calculate particle weight:
                    const double GSUM  = ion.p_BC.GSUM;
                    const double a_p_new = min(GSUM/uN_total, 1000.0);

                    ion.p_BC.a_p_new = a_p_new;

                    // Reset accumulators:
                    ion.p_BC.S1   = 0;
                    ion.p_BC.S2   = 0;
                    ion.p_BC.GSUM = 0;
                }

            } // BC

        } //  Species

    } // Particle MPIs

}

void particleBC_TYP::getFluxesAcrossBoundaries(const params_TYP &params, const CS_TYP &CS, fields_TYP &fields, vector<ionSpecies_TYP> &IONS)
{
    if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        dot_.N1 = 0;
        dot_.E1 = 0;
        dot_.N2 = 0;
        dot_.E2 = 0;

        // Simulation time step:
        const double DT = params.DT;

        for (ionSpecies_TYP &ion : IONS)
        {
            const double NCP = ion.NCP/DT;

            const int iie=ion.NSP;
            #pragma omp declare reduction(sum : struct dot_buffer : omp_out.N1 += omp_in.N1, omp_out.E1 += omp_in.E1, omp_out.N2 += omp_in.N2, omp_out.E2 += omp_in.E2)
            #pragma omp parallel for default(none) shared(params, ion, NCP, iie) reduction(sum:dot_)
            for(int ii=0; ii<iie; ii++)
            {
                // Boundary 1:
                if ( ion.f1(ii) == 1 )
                {
                    const double a_p = ion.a_p(ii);

                    // Accumulate fluxes:
                    dot_.N1 += NCP*a_p;
                    dot_.E1 += NCP*a_p*ion.dE1(ii);

                }

                // Boundary 2:
                if ( ion.f2(ii) == 1 )
                {
                    const double a_p = ion.a_p(ii);

                    // Accumulate fluxes:
                    dot_.N2 += NCP*a_p;
                    dot_.E2 += NCP*a_p*ion.dE2(ii);

                } // if

            } // omp parallel for

        } // Species

        // Reduce over all MPI process
        MPI_AllreduceDouble<4> (params,&dot_.N1);

    } // Particle MPI
}

// =============================================================================
void particleBC_TYP::applyParticleReinjection(const params_TYP &params, const CS_TYP &CS, fields_TYP &fields, vector<ionSpecies_TYP> &IONS)
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
    if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        for (ionSpecies_TYP &ion : IONS)
        {

            // Ion mass:
            const double Ma = ion.M*0.5;

            #pragma omp parallel default(none) shared(params, CS, fields, ion, std::cout, Ma)
            {
                uniform_2Pi &rand_2pi = randoms_2pi[picos::random::thread()];
                uniform_one &rand_one = randoms_one[picos::random::thread()];

                const int iie=ion.NSP;
                #pragma omp parallel for
                for(int ii=0; ii<iie; ii++)
                {
                    if ( ion.f1(ii) == 1 || ion.f2(ii) == 1 )
                    {
                        // Re-inject particle:
                        // ===================
                        particleReinjection(ii, params, CS, fields, ion, rand_2pi, rand_one);
                        
                        // Newly injected flag:
                        // ====================
                        ion.f5(ii)  = 1;
                        ion.dE5(ii) = Ma*dot(ion.V_p.row(ii), ion.V_p.row(ii));
                        
                        // Reset injection flag:
                        // =====================
                        ion.f1(ii) = 0;
                        ion.f2(ii) = 0;
                        
                        // Reset Exit energy:
                        // =====================
                        ion.dE1(ii) = 0;
                        ion.dE2(ii) = 0;
                        
                    } // flag guard
                } // pragma omp for

            } // pragma omp parallel

    	} //  Species

	} // Particle MPIs

    // Calculate actual fueling rate and power:
    // =======================================================
    getParticleInjectionRates(params,CS,fields,IONS);

}

void particleBC_TYP::getParticleInjectionRates(const params_TYP &params, const CS_TYP &CS, fields_TYP &fields, vector<ionSpecies_TYP> &IONS)
{
    if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        dot_.N5 = 0;
        dot_.E5 = 0;
        const double DT = params.DT;

        for (ionSpecies_TYP &ion : IONS)
        {
            const double NCP = ion.NCP/DT;

            const int iie=ion.NSP;
            #pragma omp declare reduction(sum : struct dot_buffer : omp_out.N5 += omp_in.N5, omp_out.E5 += omp_in.E5)
            #pragma omp parallel for default(none) shared(ion, iie, NCP) reduction(sum:dot_)
            for(int ii=0; ii<iie; ii++)
            {
                // Boundary 1:
                if (ion.f5(ii) == 1 )
                {
                    const double a_p = ion.a_p(ii)*NCP;

                    // Accumulate fluxes:
                    dot_.N5 += a_p;
                    dot_.E5 += a_p*ion.dE5(ii);

                    // Clear flags:
                    ion.f5(ii)  = 0;
                    ion.dE5(ii) = 0;
                }
            } // omp parallel for

        } // Species

        // Reduce over all MPI process
        MPI_AllreduceDouble<2> (params,&dot_.N5);

    } // Particle MPI
}

// =============================================================================
void particleBC_TYP::particleReinjection(const int ii, const params_TYP &params, const CS_TYP &CS, fields_TYP &fields, ionSpecies_TYP &ION, uniform_2Pi &rand_2pi, uniform_one &rand_one) const
{
    // Particle velocity:
	// =========================================================================

    // 1: Warm plasma source, 2: NBI or 4: simple-reinjection:
	if (ION.p_BC.BC_type == 1 || ION.p_BC.BC_type == 2 || ION.p_BC.BC_type == 4)
	{
		double T;
		double E;

		if (ION.p_BC.BC_type == 1 || ION.p_BC.BC_type == 4) // Warm plasma source
		{
			T = ION.p_BC.T;
			E = 0;

            // cout << "T = " << T*CS->temperature << endl;
		}
		if (ION.p_BC.BC_type == 2) // NBI
		{
			T = ION.p_BC.T;
			E = ION.p_BC.E;
		}

		// Mass of ion:
		const double Ma = ION.M;

		// Thermal velocity of source:
		const double vT = sqrt(2*F_E_DS*T/Ma);

		// Pitch angle of source:
		const double xip = cos(ION.p_BC.eta);

		// Drift velocity of source:
		const double U  = sqrt(2*F_E_DS*E/Ma);
		const double Ux = U*xip;
		const double Uy = U*sqrt(1 - xip*xip);
		const double Uz = 0;

		// Thermal spread: Note I have prefactored out a 1/Sqrt(2) and removed
        // Sqrt(2) term from the R_1 and R_2.
		const double sigma_v = vT;

		// Box muller:
		const double R_1 = sigma_v*sqrt(-log(rand_one()));
        const double t_2 = rand_2pi();
        const double R_3 = sigma_v*sqrt(-log(rand_one()));
        const double t_4 = rand_2pi();

		// Thermal component:
        const double wx = R_3*cos(t_4);
        const double wy = R_1*cos(t_2);
        const double wz = R_1*sin(t_2);

		// Total velocity components:
        const double v_par = Ux + wx;
        const double v_y   = Uy + wy;
        const double v_z   = Uz + wz;
        const double v_per = hypot(v_y, v_z);

        // Assign to IONS:
        ION.V_p(ii,0) = v_par;
        ION.V_p(ii,1) = v_per;

        // cout << "vpar = " << v_par*CS->velocity << endl;
        // cout << "vper = " << v_per*CS->velocity << endl;

	}

    // Periodic:
//    if (ION.p_BC.BC_type == 3 )
//    {
        // Do nothing
//    }

	// Particle position:
    // =========================================================================
    // 1: Warm plasma source, 2: NBI or 4: simple-reinjection:
	if (ION.p_BC.BC_type == 1 || ION.p_BC.BC_type == 2 || ION.p_BC.BC_type == 4)
	{
		// Gaussian distribution in space:
		const double mean_x  = ION.p_BC.mean_x;
		const double sigma_x = ION.p_BC.sigma_x*numbers::sqrt2_v<double>;
        double new_x = mean_x + sigma_x*sqrt(-log(rand_one()))*cos(rand_2pi());

        // Domain boundaries:
        const double LX_min = params.geometry.LX_min;
        const double LX_max = params.geometry.LX_max;

        // Correction to prevent injecting out of bounds:
        while( (new_x < LX_min) || (new_x > LX_max) )
		{
			 //cout<<"Out of bound, Xp(ii) = "<< new_x << endl;

             // Assign new postion:
			 new_x = mean_x + sigma_x*sqrt(-log(rand_one()))*cos(rand_2pi());

//             if ( (new_x > LX_min) || (new_x < LX_max) )
//             {
                 //cout<< "Corrected X = " << new_x <<  endl;
//             }
		}

		ION.X_p(ii) = new_x;
	}

    //4: Periodic:
	if (ION.p_BC.BC_type == 3)
	{
        if (ION.X_p(ii) > params.geometry.LX_max)
        {
            ION.X_p(ii) = params.geometry.LX_min;
        }

        if (ION.X_p(ii) < params.geometry.LX_min)
        {
            ION.X_p(ii) = params.geometry.LX_max;
        }

	}

	// Particle weight:
    // =========================================================================
    // 1: Warm plasma source or 2: NBI
	if (ION.p_BC.BC_type == 1 || ION.p_BC.BC_type == 2)
	{
		ION.a_p(ii) = ION.p_BC.a_p_new;
	}

    // 3: Periodic or 4: basic re-injection
//    if (ION.p_BC.BC_type == 3 || ION.p_BC.BC_type == 4)
//    {
//        // Do nothing
//        //IONS->a_p(ii) = 1;
//    }

}
