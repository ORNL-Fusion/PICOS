#include "initDistribution.h"

namespace
{
std::uint_fast64_t initialConditionSeed(const params_TYP *params, const int speciesIndex)
{
    return static_cast<std::uint_fast64_t>(params->initialConditionRandomSeed)
        + 104729ULL*static_cast<std::uint_fast64_t>(speciesIndex + 1)
        + 4099ULL*static_cast<std::uint_fast64_t>(params->mpi.MPI_DOMAIN_NUMBER + 1);
}
}

initDist_TYP::initDist_TYP(const params_TYP * params)
{

}


//This function creates a Maxwellian velocity distribution for IONS with a homogeneous spatial distribution.
void initDist_TYP::uniform_maxwellianDistribution(const params_TYP * params, ionSpecies_TYP * IONS, const int speciesIndex)
{
    // Uniformely distribute positions along domain:
    if (params->initialConditionRandomSeed >= 0)
    {
        arma_rng::set_seed(initialConditionSeed(params, speciesIndex));
    }
    else
    {
        arma_rng::set_seed_random();
    }
    IONS->X_p = params->geometry.LX_min + randu<vec>(IONS->NSP)*(params->geometry.LX_max - params->geometry.LX_min);

    // Maxwellian distribution for the velocity using Box-Muller:
    //arma_rng::set_seed_random();
	arma::vec R = randu(IONS->NSP);
	//arma_rng::set_seed_random();
	arma::vec phi = 2.0*M_PI*randu<vec>(IONS->NSP);

	arma::vec V2 = IONS->VTper*sqrt( -log(1.0 - R) ) % cos(phi);
	arma::vec V3 = IONS->VTper*sqrt( -log(1.0 - R) ) % sin(phi);
    if (params->velocityDistributionModel == VELOCITY_DISTRIBUTION_FORTRAN_CORRELATED_PERP)
    {
        V3 = V2;
    }
    arma::vec V4 = sqrt( pow(V2,2) + pow(V3,2) );

	//arma_rng::set_seed_random();
	R = randu<vec>(IONS->NSP);
	//arma_rng::set_seed_random();
	phi = 2.0*M_PI*randu<vec>(IONS->NSP);

	arma::vec V1 = IONS->VTpar*sqrt( -log(1.0 - R) ) % sin(phi);

    // Assign velocities:
    IONS->V_p.col(0) = V1;
    if (params->advanceParticleMethod == PARTICLE_PUSH_BORIS_FULL_ORBIT)
    {
        IONS->V_p.col(1) = V2;
        IONS->V_p.col(2) = V3;
    }
    else
    {
        IONS->V_p.col(1) = V4;
    }
}

void initDist_TYP::nonuniform_maxwellianDistribution(const params_TYP * params, ionSpecies_TYP * IONS, const int speciesIndex)
{
    if (params->initialConditionRandomSeed >= 0)
        arma_rng::set_seed(initialConditionSeed(params, speciesIndex));
    else
        arma_rng::set_seed_random();

    const arma::vec &x = IONS->p_IC.x_profile;
    arma::vec spatialPdf = IONS->p_IC.densityFraction_profile
                         % (params->em_IC.BX/params->em_IC.Bx_profile);
    arma::vec cdf(x.n_elem, fill::zeros);
    for (arma::uword ii=1; ii<x.n_elem; ++ii)
        cdf(ii) = cdf(ii-1) + 0.5*(spatialPdf(ii-1) + spatialPdf(ii))*(x(ii) - x(ii-1));
    if (!(cdf(cdf.n_elem-1) > 0.0))
        throw std::runtime_error("Nonuniform initial-condition density integral must be positive");
    cdf /= cdf(cdf.n_elem-1);

    arma::vec X;
    interp1(cdf, x, randu<vec>(IONS->NSP), X, "linear");

    arma::vec TparLocal;
    arma::vec TperLocal;
    interp1(x, IONS->p_IC.Tpar_profile, X, TparLocal, "linear");
    interp1(x, IONS->p_IC.Tper_profile, X, TperLocal, "linear");
    const arma::vec sigmaPar = (IONS->VTpar/std::sqrt(2.0))
                             * sqrt(TparLocal/IONS->p_IC.Tpar);
    const arma::vec sigmaPer = (IONS->VTper/std::sqrt(2.0))
                             * sqrt(TperLocal/IONS->p_IC.Tper);
    const arma::vec V1 = sigmaPar%randn<vec>(IONS->NSP);
    const arma::vec V2 = sigmaPer%randn<vec>(IONS->NSP);
    const arma::vec V3 = sigmaPer%randn<vec>(IONS->NSP);

    // Assign value to "V":
    // ====================
    arma::vec V4 = sqrt( pow(V2,2) + pow(V3,2) );
    IONS->V_p.col(0) = V1;
    if (params->advanceParticleMethod == PARTICLE_PUSH_BORIS_FULL_ORBIT)
    {
        IONS->V_p.col(1) = V2;
        IONS->V_p.col(2) = V3;
    }
    else
    {
        IONS->V_p.col(1) = V4;
    }

    IONS->X_p = X;

}
