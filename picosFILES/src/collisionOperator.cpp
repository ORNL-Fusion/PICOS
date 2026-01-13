#include <cmath>
#include <numbers>
#include <algorithm>

#include "collisionOperator.h"

using namespace std;

// Velocity scattering operator:
// =============================================================================
void coll_operator_TYP::u_CollisionOperator(double &w,
                                            const double xab,
                                            const double wTb,
                                            const double nb,
                                            const double Tb,
                                            const double Mb,
                                            const double Zb,
                                            const double Za,
                                            const double Ma,
                                            const double DT,
                                            uniform_random &rand)
{
    const double BoozerFactor = 1.0;
    const uint8_t energyOperatorModel = 2;

    // Normalized collision rate:
    double nu_E_dt = BoozerFactor*nu_E<energyOperatorModel> (xab,nb,Tb,Mb,Zb,Za,Ma)*DT;

    // Calculate substeps:
    const int Nstep = std::min(static_cast<int> (nu_E_dt*2.5) + 1, 100);

    // Apply operator:
    nu_E_dt  = nu_E_dt/Nstep;

    const double mof = Ma/(2*F_E);
    const double B = 2.0*nu_E_dt*( 1.5 + E_nuE_d_nu_E_dE(xab))*Tb;
    const double Tbnu_e_dt = Tb*nu_E_dt;
    const double twonu_e_dt = -2*nu_E_dt;

    w = w*w;
    for (int kk = 0; kk<Nstep; kk++)
    {
        const double E0 = mof*w;
        const double A = twonu_e_dt*E0;

        // Random number between 0 and 1:
        const short Rm = 2*rand() - 1;

        const double C = 2*Rm*sqrt(Tbnu_e_dt*E0);
        w = (E0 + A + B + C)/mof;
    }
    w = sqrt(w);
}

// Pitch angle scattering operator:
// =============================================================================
void coll_operator_TYP:: xi_CollisionOperator(double &xi,
                                              const double xab,
                                              const double wTb,
                                              const double nb,
                                              const double Tb,
                                              const double Mb,
                                              const double Zb,
                                              const double Za,
                                              const double Ma,
                                              const double DT,
                                              uniform_random &rand)
{
    // Normalized collisional rate:
    // ===========================
    double nu_D_dt = nu_D(xab,nb,Tb,Mb,Zb,Za,Ma)*DT;

// NOTE: The xi collision operator computes nu_D_dt before truncating Nstep so
//       We cannot make nstep const here but we could in the u operator.
    
    // Calculate substeps:
    // ===========================
    int Nstep   = round(nu_D_dt*2.5) + 1;

    // Recalculate normalized rate:
    // ============================
    nu_D_dt  = nu_D_dt/Nstep;

    // Limit substepping:
    // ===========================
    if (Nstep > 100)
    {
        cout << "Nstep for 'xi' operator = " << Nstep << endl;
        Nstep = 100;
    }

    // Apply operator:
    // ===========================
    for (int kk = 0; kk<Nstep; kk++)
    {
        // Deterministic part:
        // ==================
        const double A = -xi*nu_D_dt;

        // Stochastic part:
        // ===============
        // Random number between 0 and 1:
        const short Rm = 2*rand() - 1;

        const double C = Rm*sqrt((1.0 - xi*xi)*nu_D_dt);

        // Monte-Carlo change:
        // ==================
       xi += A + C;
    }

}

// Interpolate ion moments: test species "a" and background species "b"
void coll_operator_TYP::interpolateIonMoments(const params_TYP &params, ionSpecies_TYP &iona, const ionSpecies_TYP &ionb) const
{
    // Interpolate:
    interpolateScalarField(params, iona, ionb.n_m   , iona.n_p   );
    interpolateScalarField(params, iona, ionb.nv_m  , iona.nv_p  );
    interpolateScalarField(params, iona, ionb.Tpar_m, iona.Tpar_p);
    interpolateScalarField(params, iona, ionb.Tper_m, iona.Tper_p);
}

// Interpolate electron temperature : test species "a" and electron fluid as species "b"
void coll_operator_TYP::interpolateElectronTemperature(const params_TYP &params, ionSpecies_TYP &ion, const electrons_TYP &electrons) const
{
    // Interpolate:
    interpolateScalarField(params, ion, electrons.Te_m, ion.Te_p);
}

// Fill ghost cells:
// =============================================================================
void coll_operator_TYP::fill4Ghosts(arma::vec &v) const
{
	const int N = v.n_elem;

    v.subvec(N-2,N-1) = v.subvec(N-4,N-3);
    v.subvec(0,1)     = v.subvec(2,3);
}

// General scalar field second order interpolation method:
// =============================================================================
void coll_operator_TYP::interpolateScalarField(const params_TYP &params, const ionSpecies_TYP &ion, const arma::vec &F_m, arma::vec &F_p) const
{
	const int NX =  params.mesh.NX_IN_SIM + 4; //Mesh size along the X axis (considering the gosht cell)

	arma::vec F = zeros(NX);

	F.subvec(1,NX-2) = F_m;

	fill4Ghosts(F);

	#pragma omp parallel for default(none) shared(params, ion, F_p, F)
	for(int ii=0, iie = ion.NSP; ii<iie; ii++)
	{
		const int ix = ion.mn(ii) + 2;

		F_p(ii) += ion.wxl(ii)*F(ix-1);
		F_p(ii) += ion.wxc(ii)*F(ix);
		F_p(ii) += ion.wxr(ii)*F(ix+1);

	}//End of the parallel region
}

// Entire collision operator method:
// =============================================================================
void coll_operator_TYP::ApplyCollisions_AllSpecies(const params_TYP &params, const CS_TYP &CS, vector<ionSpecies_TYP> &IONS, electrons_TYP &electrons)
{
    if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        // Number of ION species:
    	// =====================
    	const int numIonSpecies = IONS.size();
        const int bbe = numIonSpecies + 1;

        // Time step:
        // =========
        const double DT = params.DT*CS.time;

    	for (ionSpecies_TYP &iona : IONS)
    	{
            // Number of particles is "aa" species:
        	// ===================================
        	const int NSP_a = iona.NSP;

        	// Species "aa" parameters:
        	// =======================
        	const double Ma = iona.M*CS.mass;
        	const double Za = iona.Z;

            // Initialize Species "bb" parameters:
            // =======================
            double Mb = 0.0;
            double Zb = 0.0;
            arma::vec nb  =  zeros(NSP_a,1);
            arma::vec Tb  =  zeros(NSP_a,1);
            arma::vec uxb =  zeros(NSP_a,1);

            // Initialize total ion density and flux density:
        	// ===========================================
        	arma::vec nUx_i = zeros(NSP_a,1);
        	arma::vec n_i   = zeros(NSP_a,1);

        	for (int bb=0; bb < bbe; bb++)
            {
                // Background species "bb" conditions:
				// ==================================
				if (bb < numIonSpecies) // Ions:
				{
					// Background parameters:
                    ionSpecies_TYP &ionb = IONS[bb];
					Mb = ionb.M*CS.mass;
					Zb = ionb.Z;

					// Interpolate moments:
					interpolateIonMoments(params, iona, ionb);
					const arma::vec nv_p   = iona.nv_p*CS.velocity/CS.volume;

					// Background conditions:
					nb = iona.n_p/CS.volume;
					Tb = 0.5*(iona.Tpar_p + iona.Tper_p)*CS.temperature*F_KB/F_E;
					uxb = nv_p/nb;

					// Accumulate total ion density and ion flux density:
					n_i   = n_i + nb*Zb;
					nUx_i = nUx_i + nv_p*Zb;
				}
				else // Electrons:
				{
					// Background parameters:
					Mb = F_ME;
					Zb = -1;

                    // Interpolate electron temperature:
                    interpolateElectronTemperature(params,iona,electrons);
                    Tb = iona.Te_p*CS.temperature*F_KB/F_E;

					// Background conditions:
					nb  = n_i;
					uxb = nUx_i/n_i;
				}

                // Apply collisions to all particles:
				// ==================================
                #pragma omp parallel default(none) shared(params, iona, CS, Ma, Za, Mb, Zb, nb, Tb, uxb, DT, std::cout, NSP_a)
                {
                    uniform_random &rand = randoms[picos::random::thread()];

                    #pragma omp for firstprivate(NSP_a)
                    for(int ii=0; ii<NSP_a; ii++)
                    {
                        // Species "aa":
                        // =============================================================================
                        // Velocities:
                        double vxa = iona.V_p(ii,0)*CS.velocity;
                        double vya = iona.V_p(ii,1)*CS.velocity;
                        double vza = 0;
                        
                        // Convert to ion species "bb" frame:
                        // =============================================================================
                        double wxa = vxa - uxb(ii);
                        double wya = vya;
                        double wza = vza;
                        
                        // Convert velocity from cartesian to spherical coordinate system:
                        // =============================================================================
                        double w;
                        double xi;
                        double phi;
                        cartesian2Spherical(wxa, wya, wza, w, xi, phi);
                        
                        // Apply Monte-Carlo collision operator:
                        // =============================================================================
                        double phi0 = phi;
                        
                        double wTb = sqrt(2*F_E*Tb(ii)/Mb);
                        double xab = w/wTb;
                        
                        // Velocity operator:
                        u_CollisionOperator(w, xab, wTb, nb(ii), Tb(ii), Mb, Zb, Za, Ma, DT, rand);

                        // Pitch angle operator:
                        xi_CollisionOperator(xi, xab, wTb, nb(ii), Tb(ii), Mb, Zb, Za, Ma, DT, rand);
                        
                        // Final Velocity:
                        // =============================================================================
                        // Final pitch angle:
                        // =============================================================================
                        // Reflective boundary condition:
                        xi = xi*xi > 1 ? copysign(1,xi) - fmod(xi, copysign(1,xi)) : xi;
                        
                        // Convert velocity from spherical to cartesian coordinate sytem:
                        // =====================================================================
                        Spherical2Cartesian(w, xi, phi, wxa, wya, wza);
                        
                        // Back to lab frame and normalize:
                        // =====================================================================
                        iona.V_p(ii,0) = (wxa + uxb(ii))/CS.velocity;
                        iona.V_p(ii,1) = wya/CS.velocity;
                        
                        if ( isnan(iona.V_p(ii,0)) || isnan(iona.V_p(ii,1)) )
                        {
                            cout << "isnan(V) == 1" << endl;
                        }
                        
                    } // "ii" particle loop
                }

            } // "bb" species loop

        } // "aa" species loop

    } // MPI if statement

}

// Coordinate transformation function:
// =============================================================================
void coll_operator_TYP::cartesian2Spherical(const double wx, const double wy, const double wz, double &w, double &xi, double &phi) const
{
    w = hypot(wx, wy, wz);
    xi = wx/w;
    phi = atan2(-wy,wz);
}

void coll_operator_TYP::Spherical2Cartesian(const double w, const double xi, const double phi, double &wx, double &wy, double &wz) const
{
    const double wper = w*sqrt(1.0 - xi*xi);
    wx   = w*xi;
    wy   = -wper*sin(phi);
    wz   = +wper*cos(phi);
}

double coll_operator_TYP::nu_D(const double xab, const double nb, const double Tb, const double Mb, const double Zb, const double Za, const double Ma) const
{
    return nu_ab0(nb,Tb,Mb,Zb,Za,Ma)*(erf(xab) - Gb(xab))/(xab*xab*xab);
}

double coll_operator_TYP::nu_ab0(const double nb, const double Tb, const double Mb, const double Zb, const double Za, const double Ma) const
{
    const double wTb = sqrt(2.0*F_E*Tb/Mb);
    const double wTb3 = wTb*wTb*wTb;
    const double F_E4 = F_E*F_E*F_E*F_E;
    const double ZaZb2 = Za*Zb*Za*Zb;
    return nb*F_E4*ZaZb2*logA(nb,Tb)/(2.0*numbers::pi_v<double>*Ma*Ma*F_EPSILON*F_EPSILON*wTb3);
}

double coll_operator_TYP::logA(const double nb, const double Tb) const
{
    const double Tb_sr = sqrt(Tb);
    return 30.0 - 0.5*log(nb/(Tb_sr*Tb_sr*Tb_sr));
}

double coll_operator_TYP::Gb(const double xab) const
{
    if (xab < 0.01)
    {
        return (2.0*numbers::inv_sqrtpi_v<double>/3)*xab;
    }
    else
    {
        return (erf(xab) - xab*erfp(xab))/(2.0*xab*xab);
    }
}

double coll_operator_TYP::erfp(const double xab) const
{
    return 2.0*numbers::inv_sqrtpi_v<double>*exp(-xab*xab);
}

double coll_operator_TYP::erfpp(double xab) const
{
    return -2.0*xab*erfp(xab);
}

double coll_operator_TYP::E_nuE_d_nu_E_dE(const double xab) const
{
    return 0.5*((3.0*(xab*erfp(xab) - erf(xab)) - xab*xab*erfpp(xab))/(erf(xab) - xab*erfp(xab)));
}
