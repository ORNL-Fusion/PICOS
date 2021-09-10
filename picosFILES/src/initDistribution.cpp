#include "initDistribution.h"

initDist_TYP::initDist_TYP(const params_TYP * params)
{

}


//This function creates a Maxwellian velocity distribution for IONS with a homogeneous spatial distribution.
void initDist_TYP::uniform_maxwellianDistribution(const params_TYP * params, ionSpecies_TYP * IONS)
{
    // Uniformely distribute positions along domain:
    arma_rng::set_seed_random();
    IONS->X_p = params->geometry.LX_min + randu<vec>(IONS->NSP)*(params->geometry.LX_max - params->geometry.LX_min);

    // Maxwellian distribution for the velocity using Box-Muller:
    //arma_rng::set_seed_random();
	arma::vec R = randu(IONS->NSP);
	//arma_rng::set_seed_random();
	arma::vec phi = 2.0*M_PI*randu<vec>(IONS->NSP);

	arma::vec V2 = IONS->VTper*sqrt( -log(1.0 - R) ) % cos(phi);
	arma::vec V3 = IONS->VTper*sqrt( -log(1.0 - R) ) % sin(phi);
    arma::vec V4 = sqrt( pow(V2,2) + pow(V3,2) );

	//arma_rng::set_seed_random();
	R = randu<vec>(IONS->NSP);
	//arma_rng::set_seed_random();
	phi = 2.0*M_PI*randu<vec>(IONS->NSP);

	arma::vec V1 = IONS->VTpar*sqrt( -log(1.0 - R) ) % sin(phi);

    // Assign velocities:
    IONS->V_p.col(0) = V1;
    IONS->V_p.col(1) = V4;
}

double initDist_TYP::target(const params_TYP * params,  ionSpecies_TYP * IONS, double X, double V3, double V2, double V1)
{
    // Profile for interpolation:
    arma::vec x_profile    = IONS->p_IC.x_profile;
    arma::vec Tpar_profile = IONS->p_IC.Tpar_profile;
    arma::vec Tper_profile = IONS->p_IC.Tper_profile;
    arma::vec densityFraction_profile = IONS->p_IC.densityFraction_profile;

    // Uncompress the density:
    arma::vec Bx_profile = params->em_IC.Bx_profile;
    double Bx0 = params->em_IC.BX;
    densityFraction_profile = densityFraction_profile%pow(Bx0/Bx_profile,1.5);

    /*
    if (params->mpi.IS_PARTICLES_ROOT)
    {
        cout << "################" <<  endl;
        cout << "################" <<  endl;
        cout << densityFraction_profile <<  endl;
        cout << "################" <<  endl;
        cout << "################" <<  endl;
    }
    */

    // Axial position query point:
    arma::vec xx(1,1);
    xx(0,0)= X;

    arma::vec TT3(1,1);
    interp1(x_profile,Tpar_profile,xx,TT3);
    double T3 = TT3(0,0);   //Temperature profile in x

    arma::vec TT2(1,1);
    interp1(x_profile,Tper_profile,xx,TT2);
    double T2 = TT2(0,0);   //Temperature profile in y

    arma::vec TT1(1,1);
    interp1(x_profile,Tper_profile,xx,TT1);
    double T1 = TT1(0,0);   //Temperature profile in z

    double k3=sqrt((IONS->M)/(2.0*M_PI*F_KB*T3));
    double k2=sqrt((IONS->M)/(2.0*M_PI*F_KB*T2));
    double k1=sqrt((IONS->M)/(2.0*M_PI*F_KB*T1));
    double s3=sqrt((IONS->M)/(2.0*F_KB*T3));
    double s2=sqrt((IONS->M)/(2.0*F_KB*T2));
    double s1=sqrt((IONS->M)/(2.0*F_KB*T1));

    double h= (k1*k2*k3)*exp(-sqrt(2)*(((V1*V1)*(s1*s1))+((V2*V2)*(s2*s2))+((V3*V3)*(s3*s3)))); // Pdf for 3-Velocities with temperature profile

    arma::vec gg(1,1);
    interp1(x_profile,densityFraction_profile,xx,gg); //Ne is multiplied with the compression factor
    double g = gg(0,0); //density profile

    return(g*h); //target 4-D Pdf
}


void initDist_TYP::nonuniform_maxwellianDistribution(const params_TYP * params, ionSpecies_TYP * IONS)
{
    // Initialize ion variables:
    // =========================
    //IONS->X_p = zeros(IONS->NSP);
    //IONS->V_p = zeros(IONS->NSP,2);

    // Apply MH algorith for sampling the 4-D target PDF:
    // ==================================================
    arma::vec X(IONS->NSP,fill::zeros);  // Particles position along x-axis that will be generated using M-H.
    arma::vec V3(IONS->NSP,fill::zeros); //Velocity profile in X
    arma::vec V2(IONS->NSP,fill::zeros); //Velocity profile in Y
    arma::vec V1(IONS->NSP,fill::zeros); //Velocity profile in Z

    // Initial state of search:
    // ========================
    double X_test = 0.0;
    double V3_test= 0.0;
    double V2_test= 0.0;
    double V1_test= 0.0;

    // Seed the random number generator:
    // ================================
    std::default_random_engine generator( (params->mpi.MPI_DOMAIN_NUMBER + 1 + time(NULL))*1000 );

    // Create uniform random number generator in [0,1]:
    // ================================================
    std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

    // Test number number generator:
    // ============================
    bool flagTest = false;
    if (flagTest)
    {
        std::cout << "From Rank: " << params->mpi.MPI_DOMAIN_NUMBER << std::endl;
        std::cout << uniform_distribution(generator) << " " << std::endl;
    }

    // Apply MH algorithm:
    // ====================
    unsigned int iterator = 1;
    double ratio = 0.0;
    X(0)  = params->geometry.LX_min + 0.5*(params->geometry.LX_max - params->geometry.LX_min);
    V3(0) = 0.25*IONS->VTpar;
    V2(0) = 0.1*IONS->VTpar;
    V1(0) = 0.5*IONS->VTpar;

    // Iterate over all particles in process:
    while(iterator < IONS->NSP)
    {
        // Select search point in phase space:
        X_test  = params->geometry.LX_min + uniform_distribution(generator)*(params->geometry.LX_max - params->geometry.LX_min);
        V3_test = 10*(IONS->VTper)*((uniform_distribution(generator))-0.5);
        V2_test = 10*(IONS->VTper)*((uniform_distribution(generator))-0.5);
        V1_test = 10*(IONS->VTpar)*((uniform_distribution(generator))-0.5);

        // Calculate the probability ratio:
        ratio = target(params, IONS, X_test,V3_test, V2_test, V1_test)/target(params, IONS,  X(iterator - 1),V3(iterator - 1),V2(iterator - 1),V1(iterator - 1));

        // Choose outcome based on "ratio":
        if (ratio >= 1.0)
        {
            X(iterator) = X_test;
            V3(iterator) = V3_test;
            V2(iterator) = V2_test;
            V1(iterator) = V1_test;
            iterator += 1;
        }
        else
        {
          if (uniform_distribution(generator) < ratio)
          {
              X(iterator) = X_test;
              V3(iterator) = V3_test;
              V2(iterator) = V2_test;
              V1(iterator) = V1_test;
              iterator += 1;
          }
        }
    }

    // Assign value to "V":
    // ====================
    arma::vec V4 = sqrt( pow(V2,2) + pow(V3,2) );
    IONS->V_p.col(0) = V1;
    IONS->V_p.col(1) = V4;

    // Assign value to "x":
    // ====================
    IONS->X_p = X;

}
