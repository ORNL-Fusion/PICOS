#include <cmath>
#include <numbers>
#include <algorithm>
#include <cassert>

#include "collisionOperator.h"

using namespace std;

// Velocity scattering operator:
// =============================================================================
void coll_operator_TYP::u_CollisionOperator(double &w,
                                            const double xab,
                                            const double xab2,
                                            const double nuab0,
                                            const double Tb,
                                            const double Mb,
                                            const double Ma,
                                            const double erf_xab,
                                            const double erfp_xab,
                                            const double xerfp_xab,
                                            const double gb,
                                            const double DT,
                                            const int collOperType,
                                            uniform_random &randuni)
{
    const uint8_t energyOperatorModel = 2;
    const double massTerm = (collOperType == COLLISION_OPERATOR_BOOZER_KIM) ? (1.0 + Mb/Ma) : 1.0;

    // Normalized collision rate:
    double nu_E_dt = nu_E<energyOperatorModel> (xab,nuab0,erfp_xab,Mb,Ma,gb)*DT/massTerm;
    if (!isfinite(nu_E_dt) || nu_E_dt <= 0.0)
    {
        return;
    }

    // Calculate substeps:
    const double uSubstepsRaw = min(max(nu_E_dt*2.5, 0.0), 99.0);
    const size_t Nstep = static_cast<size_t>(round(uSubstepsRaw)) + 1;

    // Apply operator:
    nu_E_dt = nu_E_dt/Nstep;

    const double mof = Ma/(2*F_E);
    const double B = 2.0*nu_E_dt*( 1.5 + E_nuE_d_nu_E_dE(xab2, erf_xab, xerfp_xab))*Tb;
    const double Tbnu_e_dt = Tb*nu_E_dt;
    const double A = 1 - 2*nu_E_dt;

    w = max(w*w, 0.0);
    for (size_t kk = 0; kk<Nstep; kk++)
    {
        const double E0 = max(mof*w, 0.0);

        // Random number ±2:
        const short Rm = 4*randuni() - 2;

        const double C = Rm*sqrt(max(0.0, Tbnu_e_dt*E0));
        // w is actually w^2 here. Clamp negative stochastic excursions to zero.
        w = max((E0*A + B + C)/mof, 0.0);
    }
    w = sqrt(w);
}

// Pitch angle scattering operator:
// =============================================================================
void coll_operator_TYP:: xi_CollisionOperator(double &xi,
                                              const double xab,
                                              const double xab2,
                                              const double nuab0,
                                              const double erf_xab,
                                              const double gb,
                                              const double DT,
                                              uniform_random &randuni)
{
    // Normalized collisional rate:
    // ===========================
    double nu_D_dt = nu_D(xab, xab2, nuab0, erf_xab, gb)*DT;
    if (!isfinite(nu_D_dt) || nu_D_dt <= 0.0)
    {
        return;
    }

// NOTE: The xi collision operator computes nu_D_dt before truncating Nstep so
//       We cannot make nstep const here but we could in the u operator.

    // Calculate substeps:
    // ===========================
    const double xiSubstepsRaw = min(max(nu_D_dt*2.5, 0.0), 99.0);
    const size_t Nstep = static_cast<size_t>(round(xiSubstepsRaw)) + 1;

    // Recalculate normalized rate:
    // ============================
    nu_D_dt  = nu_D_dt/Nstep;

    // Apply operator:
    // ===========================
    for (size_t kk = 0; kk<Nstep; kk++)
    {
        if (!isfinite(xi))
        {
            xi = 0.0;
        }
        else if (xi*xi > 1.0)
        {
            xi = copysign(1.0,xi) - fmod(xi, copysign(1.0,xi));
        }

        // Deterministic part:
        // ==================
        const double A = -xi*nu_D_dt;

        // Stochastic part:
        // ===============
        // Random number between 0 and 1:
        const short Rm = 2*randuni() - 1;

        const double C = Rm*sqrt(max(0.0, 1.0 - xi*xi)*nu_D_dt);

        // Monte-Carlo change:
        // ==================
        xi += A + C;

        if (xi*xi > 1.0)
        {
            xi = copysign(1.0,xi) - fmod(xi, copysign(1.0,xi));
        }
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

        const size_t iie = ion.NSP;
	#pragma omp parallel for default(shared)
	for(size_t ii = 0; ii<iie; ii++)
	{
		const size_t ix = ion.mn(ii) + 2;

		F_p(ii)  = ion.wxl(ii)*F(ix-1);
		F_p(ii) += ion.wxc(ii)*F(ix);
		F_p(ii) += ion.wxr(ii)*F(ix+1);

}//End of the parallel region
}

void coll_operator_TYP::collisionTotals(const params_TYP &params, const CS_TYP &CS,
                                        const vector<ionSpecies_TYP> &IONS,
                                        double totals[3]) const
{
    double local[3] = {0.0, 0.0, 0.0}; // mass, parallel momentum, kinetic energy

    for (const ionSpecies_TYP &ion : IONS)
    {
        const double mass = ion.M*CS.mass;
        const int velocityColumns = ion.V_p.n_cols;
        for (size_t ii=0; ii<ion.NSP; ii++)
        {
            const double realWeight = ion.NCP*ion.a_p(ii);
            const double vx = ion.V_p(ii,0)*CS.velocity;
            double speed2 = vx*vx;
            for (int jj=1; jj<velocityColumns; jj++)
            {
                const double vv = ion.V_p(ii,jj)*CS.velocity;
                speed2 += vv*vv;
            }
            local[0] += mass*realWeight;
            local[1] += mass*realWeight*vx;
            local[2] += 0.5*mass*realWeight*speed2;
        }
    }

    MPI_Allreduce(local, totals, 3, MPI_DOUBLE, MPI_SUM, params.mpi.COMM);
}

void coll_operator_TYP::applyCollisionConservationProjection(const params_TYP &params,
                                                             const CS_TYP &CS,
                                                             vector<ionSpecies_TYP> &IONS,
                                                             const double initialTotals[3]) const
{
    double currentTotals[3] = {0.0, 0.0, 0.0};
    collisionTotals(params, CS, IONS, currentTotals);

    const double targetMass = initialTotals[0];
    const double targetMomentum = initialTotals[1];
    const double targetEnergy = initialTotals[2];
    const double currentMass = currentTotals[0];
    if (!isfinite(targetMass) || !isfinite(targetMomentum) || !isfinite(targetEnergy) ||
        !isfinite(currentMass) || targetMass <= 0.0 || currentMass <= 0.0 || targetEnergy <= 0.0)
    {
        return;
    }

    const double targetU = targetMomentum/targetMass;
    const double currentU = currentTotals[1]/currentMass;
    const double velocityShift = (targetU - currentU)/CS.velocity;
    for (ionSpecies_TYP &ion : IONS)
    {
        for (size_t ii=0; ii<ion.NSP; ii++)
        {
            ion.V_p(ii,0) += velocityShift;
        }
    }

    double shiftedTotals[3] = {0.0, 0.0, 0.0};
    collisionTotals(params, CS, IONS, shiftedTotals);
    const double targetBulkEnergy = 0.5*targetMass*targetU*targetU;
    const double targetThermalEnergy = targetEnergy - targetBulkEnergy;
    const double shiftedThermalEnergy = shiftedTotals[2] - targetBulkEnergy;
    if (!isfinite(targetThermalEnergy) || !isfinite(shiftedThermalEnergy) ||
        targetThermalEnergy <= 0.0 || shiftedThermalEnergy <= 0.0)
    {
        return;
    }

    const double scale = sqrt(targetThermalEnergy/shiftedThermalEnergy);
    if (!isfinite(scale) || scale <= 0.0)
    {
        return;
    }

    const double targetUDimensionless = targetU/CS.velocity;
    for (ionSpecies_TYP &ion : IONS)
    {
        const int velocityColumns = ion.V_p.n_cols;
        for (size_t ii=0; ii<ion.NSP; ii++)
        {
            ion.V_p(ii,0) = targetUDimensionless + scale*(ion.V_p(ii,0) - targetUDimensionless);
            for (int jj=1; jj<velocityColumns; jj++)
            {
                ion.V_p(ii,jj) *= scale;
            }
        }
    }
}

// Entire collision operator method:
// =============================================================================
void coll_operator_TYP::ApplyCollisions_AllSpecies(const params_TYP &params, const CS_TYP &CS, vector<ionSpecies_TYP> &IONS, const electrons_TYP &electrons)
{
    if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        const double tnorm = CS.temperature*F_KB/F_E;
        // Number of ION species:
        // =====================
        const size_t numIonSpecies = IONS.size();
        bool hasKineticElectrons = false;
        for (const ionSpecies_TYP &ion : IONS)
        {
            hasKineticElectrons = hasKineticElectrons || (ion.Z < 0.0);
        }
        const size_t bbe = hasKineticElectrons ? numIonSpecies : numIonSpecies + 1;

        // Time step:
        // =========
        const double DT = params.DT*CS.time;
        const bool conserveCollisionMoments = (params.SW.collisionConservationProjection == 1);
        double initialCollisionTotals[3] = {0.0, 0.0, 0.0};
        if (conserveCollisionMoments)
        {
            collisionTotals(params, CS, IONS, initialCollisionTotals);
        }

        for (ionSpecies_TYP &iona : IONS)
        {
            // Number of particles is "aa" species:
            // ===================================
            const size_t NSP_a = iona.NSP;

            // Species "aa" parameters:
            // =======================
            const double Ma = iona.M*CS.mass;
            const double Za2 = iona.Z*iona.Z;
            const bool fullOrbit = (params.advanceParticleMethod == PARTICLE_PUSH_BORIS_FULL_ORBIT && iona.V_p.n_cols > 2);

            // Initialize Species "bb" parameters:
            // =======================
            double Mb;
            double Zb2;
            arma::vec nb  =  zeros(NSP_a,1);
            arma::vec Tb  =  zeros(NSP_a,1);
            arma::vec uxb =  zeros(NSP_a,1);

            // Initialize total ion density and flux density:
            // ===========================================
            arma::vec nUx_i = zeros(NSP_a,1);
            arma::vec n_i   = zeros(NSP_a,1);

            for (size_t bb=0; bb < bbe; bb++)
            {
                // Background species "bb" conditions:
				// ==================================
				if (bb < numIonSpecies) // Ions:
				{
					// Background parameters:
                    ionSpecies_TYP &ionb = IONS[bb];
					Mb = ionb.M*CS.mass;
					Zb2 = ionb.Z*ionb.Z;

					// Interpolate moments:
					interpolateIonMoments(params, iona, ionb);
					const arma::vec nv_p   = iona.nv_p*CS.velocity/CS.volume;

					// Background conditions:
					nb = iona.n_p/CS.volume;
					Tb = 0.5*(iona.Tpar_p + iona.Tper_p)*tnorm;
					arma::vec nbSafe = nb;
					nbSafe.transform( [](double val) { return max(val, double_zero); } );
					uxb = nv_p/nbSafe;

					// Accumulate total ion density and ion flux density:
					n_i   = n_i + nb*ionb.Z;
					nUx_i = nUx_i + nv_p*ionb.Z;
				}
				else // Electrons:
				{
					// Background parameters:
					Mb = F_ME;
					Zb2 = 1;

                    // Interpolate electron temperature:
                    interpolateElectronTemperature(params,iona,electrons);
                    Tb = iona.Te_p*tnorm;

					// Background conditions:
					nb  = n_i;
					arma::vec n_i_safe = n_i;
					n_i_safe.transform( [](double val) { return max(val, double_zero); } );
					uxb = nUx_i/n_i_safe;
				}

                const double ZaZb2 = Za2*Zb2;

                // Apply collisions to all particles:
				// ==================================
                #pragma omp parallel default(none) shared(params, iona, CS, Ma, ZaZb2, Mb, nb, Tb, uxb, DT, std::cout, NSP_a) firstprivate(fullOrbit)
                {
                    uniform_random &randuni = randoms[picos::random::thread()];

                    #pragma omp for firstprivate(NSP_a)
                    for(size_t ii=0; ii<NSP_a; ii++)
                    {
                        // Species "aa":
                        // =============================================================================
                        // Velocities:
                        // Convert to ion species "bb" frame:
                        // =============================================================================
                        const double local_uxb = uxb(ii);
                        if (!isfinite(local_uxb) || !isfinite(iona.V_p(ii,0)) || !isfinite(iona.V_p(ii,1)))
                        {
                            continue;
                        }
                        if (fullOrbit && !isfinite(iona.V_p(ii,2)))
                        {
                            continue;
                        }
                        double wxa = iona.V_p(ii,0)*CS.velocity - local_uxb;
                        const double oldVy = fullOrbit ? iona.V_p(ii,1) : 0.0;
                        const double oldVz = fullOrbit ? iona.V_p(ii,2) : 0.0;
                        const double oldVper = fullOrbit ? hypot(oldVy, oldVz) : iona.V_p(ii,1);
                        double wya = oldVper*CS.velocity;

                        // Convert velocity from cartesian to spherical coordinate system:
                        // =============================================================================
                        double w;
                        double xi;
                        double sinphi;
                        cartesian2Spherical(wxa, wya, w, xi, sinphi);

                        // Apply Monte-Carlo collision operator:
                        // =============================================================================
                        const double local_tb = Tb(ii);
                        const double local_nb = nb(ii);
                        if (!isfinite(local_nb) || !isfinite(local_tb) || !isfinite(w) ||
                            local_nb <= double_zero || local_tb <= double_zero || w <= double_zero)
                        {
                            continue;
                        }
                        const double wTb = sqrt(2*F_E*local_tb/Mb);
                        if (!isfinite(wTb) || wTb <= double_zero)
                        {
                            continue;
                        }
                        const double xab = max(w/wTb, 1.0e-8);

                        // These get called multiple times in the collision
                        // operators so compute them once and pass them into
                        // the operators.
                        const double erf_xab = erf(xab);
                        // Every call to erfp has xab multiplied by it.
                        const double xab2 = xab*xab;
                        const double erfp_xab = erfp(xab2);
                        const double xerfp_xab = xab*erfp_xab;
                        const double gb = Gb(xab, xab2, erf_xab, xerfp_xab);
                        const double nuab0 = nu_ab0(wTb,nb(ii),local_tb,ZaZb2,Ma);
                        if (!isfinite(nuab0) || nuab0 <= 0.0)
                        {
                            continue;
                        }

                        // Velocity operator:
                        u_CollisionOperator(w, xab, xab2, nuab0, local_tb, Mb, Ma, erf_xab, erfp_xab, xerfp_xab, gb, DT, params.collOperType, randuni);

                        // Pitch angle operator:
                        xi_CollisionOperator(xi, xab, xab2, nuab0, erf_xab, gb, DT, randuni);

                        // Final Velocity:
                        // =============================================================================
                        // Final pitch angle:
                        // =============================================================================
                        // Reflective boundary condition:
                        xi = xi*xi > 1 ? copysign(1,xi) - fmod(xi, copysign(1,xi)) : xi;

                        // Convert velocity from spherical to cartesian coordinate sytem:
                        // =====================================================================
                        Spherical2Cartesian(w, xi, sinphi, wxa, wya);

                        // Back to lab frame and normalize:
                        // =====================================================================
                        iona.V_p(ii,0) = (wxa + local_uxb)/CS.velocity;
                        if (fullOrbit)
                        {
                            const double newVper = abs(wya/CS.velocity);
                            if (oldVper > double_zero)
                            {
                                const double scale = newVper/oldVper;
                                iona.V_p(ii,1) = oldVy*scale;
                                iona.V_p(ii,2) = oldVz*scale;
                            }
                            else
                            {
                                iona.V_p(ii,1) = newVper;
                                iona.V_p(ii,2) = 0.0;
                            }
                        }
                        else
                        {
                            iona.V_p(ii,1) = wya/CS.velocity;
                        }

                    } // "ii" particle loop
                }

            } // "bb" species loop

        } // "aa" species loop

        if (conserveCollisionMoments)
        {
            applyCollisionConservationProjection(params, CS, IONS, initialCollisionTotals);
        }

    } // MPI if statement

}

// Coordinate transformation function:
// =============================================================================
void coll_operator_TYP::cartesian2Spherical(const double wx, const double wy, double &w, double &xi, double &sinphi) const
{
    w = hypot(wx, wy);
    if (w <= double_zero)
    {
        xi = 0.0;
        sinphi = -1.0;
        return;
    }
    xi = wx/w;

//  NOTE: This is either -Pi, indeterminate, Pi depending on the value of wz so we can't
//        eliminate phi even though wz is always zero. But we can store sin(phi)
//        here istead and eliminate the calls to atan2 and sin.
//    phi = atan2(-wy, 0);
    sinphi = -copysign(1, wy);
}

void coll_operator_TYP::Spherical2Cartesian(const double w, const double xi, const double sinphi, double &wx, double &wy) const
{
    const double xiBounded = max(-1.0, min(1.0, xi));
    const double wper = w*sqrt(max(0.0, 1.0 - xiBounded*xiBounded));
    wx   = w*xiBounded;
//  NOTE: In cartesian2Spherical we eliminated the atan2 and computed sin(phi)
//        directly. So phi here is realy sign phi.
//    wy   = -wper*sin(phi);
    wy = -wper*sinphi;
}

double coll_operator_TYP::nu_D(const double xab, const double xab2, const double nuab0, const double erf_xab, const double gb) const
{
    return nuab0*(erf_xab - gb)/(xab*xab2);
}

double coll_operator_TYP::nu_ab0(const double wtb, const double nb, const double Tb, const double ZaZb2, const double Ma) const
{
    if (!isfinite(wtb) || !isfinite(nb) || !isfinite(Tb) || !isfinite(Ma) ||
        wtb <= double_zero || nb <= double_zero || Tb <= double_zero || Ma <= 0.0)
    {
        return 0.0;
    }
    const double wTb3 = wtb*wtb*wtb;
    const double F_E4 = F_E*F_E*F_E*F_E;
    return nb*F_E4*ZaZb2*logA(nb,Tb)/(2.0*numbers::pi_v<double>*Ma*Ma*F_EPSILON*F_EPSILON*wTb3);
}

double coll_operator_TYP::logA(const double nb, const double Tb) const
{
    const double nbSafe = max(nb, double_zero);
    const double TbSafe = max(Tb, double_zero);
    const double Tb_sr = sqrt(TbSafe);
    return 30.0 - 0.5*log(nbSafe/(Tb_sr*Tb_sr*Tb_sr));
}

double coll_operator_TYP::Gb(const double xab, const double xab2, const double erf_xab, const double xerfp_xab) const
{
    if (xab < 0.01)
    {
        return (2.0*numbers::inv_sqrtpi_v<double>/3)*xab;
    }
    else
    {
        return (erf_xab - xerfp_xab)/(2.0*xab2);
    }
}

double coll_operator_TYP::erfp(const double xab2) const
{
    return 2.0*numbers::inv_sqrtpi_v<double>*exp(-xab2);
}

double coll_operator_TYP::erfpp(const double xerfp_xab) const
{
    return -2.0*xerfp_xab;
}

double coll_operator_TYP::E_nuE_d_nu_E_dE(const double xab2, const double erf_xab, const double xerfp_xab) const
{
    const double denominator = erf_xab - xerfp_xab;
    if (xab2 < 1.0e-4 || abs(denominator) <= double_zero)
    {
        return 0.0;
    }
    return 0.5*((3.0*(xerfp_xab - erf_xab) - xab2*erfpp(xerfp_xab))/denominator);
}
