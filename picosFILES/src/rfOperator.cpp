#include "rfOperator.h"

RF_Operator_TYP::RF_Operator_TYP(params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {

    } // PARTICLE MPI
}

void RF_Operator_TYP::MPI_AllreduceDouble(params_TYP * params, double * v)
{
    double recvbuf = 0;

    MPI_Allreduce(v, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, params->mpi.COMM);

    *v = recvbuf;
}

void RF_Operator_TYP::calculateResNum(int ii, params_TYP * params, CS_TYP * CS, fields_TYP * fields, ionSpecies_TYP * IONS)
{
    // Ion parameters:
    double Ma = IONS->M;
    double Q  = IONS->Q;

    // Particle states:
    double vpar = IONS->V_p(ii,0);
    double Bp   = IONS->BX_p(ii);
    double wcp   = abs(Q)*Bp/Ma;

    // RF paramters:
    int n       = params->RF.n_harmonic;
    double kpar = params->RF.kpar;
    double f    = params->RF.freq;
    double wrf  = 2*M_PI*f;

    // Cyclotron resonance number:
    IONS->resNum(ii) = wrf - kpar*vpar - n*wcp;
}

void RF_Operator_TYP::calculateResNum_AllSpecies(params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
        for (int ss=0; ss<IONS->size();ss++)
        {
            int NSP   = IONS->at(ss).NSP;
            double Ma = IONS->at(ss).M;

            #pragma omp parallel for default(none) shared(params, IONS, ss, CS, fields, std::cout) firstprivate(NSP,Ma)
            for(int ii=0; ii<NSP; ii++)
            {
                // Store previous resNum:
                IONS->at(ss).resNum_(ii) = IONS->at(ss).resNum(ii);

                // Calculate new resNum:
                calculateResNum(ii,params,CS,fields,&IONS->at(ss));

            } // particles
        } // Species
}

void RF_Operator_TYP::checkResNumAndFlag_AllSpecies(params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    for (int ss=0; ss<IONS->size();ss++)
    {
        int NSP   = IONS->at(ss).NSP;

        #pragma omp parallel for default(none) shared(params, IONS, ss, CS, fields, std::cout) firstprivate(NSP)
        for(int ii=0; ii<NSP; ii++)
        {
            // Parameters needed to check resonance condition:
            double xp = IONS->at(ss).X_p(ii);
            double x1 = params->RF.x1;
            double x2 = params->RF.x2;
            double resNum  = IONS->at(ss).resNum(ii);
            double resNum_ = IONS->at(ss).resNum_(ii);

            // Check resonance condition:
            // true: argument negative
            // false: otherwise (positive or zero)
            bool dresNum_sign = signbit(resNum*resNum_);

            // Flag particles in resonance:
            if ( dresNum_sign & (xp > x1) & (xp < x2) )
            {
                IONS->at(ss).f3(ii) = 1;

                /*
                if (params->mpi.MPI_DOMAIN_NUMBER == 2)
                {
                    cout << "resNum_: " << resNum_ << endl;
                    cout << "resNum: " << resNum << endl;
                    cout << "x1: " << x1*CS->length << "[m]" << endl;
                    cout << "x2: " << x2*CS->length << "[m]" << endl;
                    cout << "xp: " << xp*CS->length << "[m]" << endl;
                    cout << "Bp: " << IONS->at(ss).BX_p(ii)*CS->bField << "[T]" << endl;
                }
                */

            }

        } // particles
    } // Species
}

void RF_Operator_TYP::calculateRfTerms(int ii, params_TYP * params, CS_TYP * CS, fields_TYP * fields, ionSpecies_TYP * IONS)
{
    // Ion parameters:
    double Ma = IONS->M;
    double Q  = IONS->Q;
    double e  = F_E_DS;
    int Z     = IONS->Z;

    // RF paramters:
    int n         = params->RF.n_harmonic;
    double kper   = params->RF.kper;
    double kpar   = params->RF.kpar;
    double tau_rf = 0;
    double rL     = 0;
    double flr    = 0;
    double J_nm1  = 0;
    double J_np1  = 0;
    double mean_dKE_per = 0;

    // Particle states:
    double vpar = IONS->V_p(ii,0);
    double vper = IONS->V_p(ii,1);

    // Particle-defined fields:
    double Bp   = IONS->BX_p(ii);
    double dBp  = IONS->dBX_p(ii);
    double ddBp = IONS->ddBX_p(ii);
    double Ep   = IONS->EX_p(ii);

    // Derived quantities:
    double KE_par  = 0.5*Ma*vper*vper;
    double KE_per  = 0.5*Ma*vper*vper;
    double Omega   = abs(Q)*Bp/Ma;
    double dOmega  = abs(Q)*dBp/Ma;
    double ddOmega = abs(Q)*ddBp/Ma;

    // Calculate the first and second time derivative of Omega:
    double Omega_dot  = vpar*dOmega;
    double Omega_ddot = pow(vpar,2)*ddOmega  - pow(vper,2)*pow(dOmega,2)/(2.*Omega)  +  (Q/Ma)*Ep*dOmega;

    // Calculate the interaction time:
    if ( pow(n*Omega_ddot,2) > 4.8175*abs(pow(n*Omega_dot,3)) )
    {
        // tau_b
        // Approximate Ai(x) ~ 0.3833
        tau_rf = (2*M_PI)*pow(abs(2/(n*Omega_ddot)),1.0/3.0)*0.3833;
        //cout << "tau_b = " << tau_rf*CS->time << "[s]" << endl;
    }
    else
    {
        // tau_a
        tau_rf = sqrt(2*M_PI/abs(n*Omega_dot));
        //cout << "tau_a = " << tau_rf*CS->time << "[s]" << endl;
    }

    // Calculate bessel terms:
    rL  = vper/Omega;
    flr = kper*rL;
    J_nm1 = std::cyl_bessel_j(n - 1,flr);
    J_np1 = std::cyl_bessel_j(n + 1,flr);

    /*
    cout << "rL = " << rL*CS->length << endl;
    cout << "FLR = " << flr << endl;
    cout << "Omega = " << Omega/CS->time << endl;
    cout << "vper = " << vper*CS->velocity << endl;
    cout << "kper = " << kper/CS->length << endl;
    cout << "J_nm1 = " << J_nm1 << endl;
    cout << "J_np1 = " << J_np1 << endl;
    */

    // Calculate the mean RF kick per unit electric field squared:
    // ********************
    if (Z > 0) // Positive ions
    {
        double E_m = 0;
        double E_p = 1;
        mean_dKE_per = 0.5*(pow(e,2)/Ma)*pow(abs(E_p*J_nm1 + E_m*J_np1)*tau_rf,2); // [J] normalized energy
    }
    if (Z < 0) // Negative particles
    {
        double E_m = 1;
        double E_p = 0;
        mean_dKE_per = 0.5*(pow(e,2)/Ma)*pow(abs(E_m*J_nm1 + E_p*J_np1)*tau_rf,2); // [J] normalized energy
    }

    // Populate output:
    IONS->udErf(ii)   = mean_dKE_per; // [J] normalized energy
    IONS->doppler(ii) = kpar*vpar/(n*Omega);
    IONS->udE3(ii)    = IONS->udErf(ii)*(1 + IONS->doppler(ii)); // [J] normalized energy
}

void RF_Operator_TYP::calculateRfTerms_AllSpecies(params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    for (int ss=0; ss<IONS->size();ss++)
    {
        int NSP   = IONS->at(ss).NSP;

        #pragma omp parallel for default(none) shared(params, IONS, ss, CS, fields, std::cout) firstprivate(NSP)
        for(int ii=0; ii<NSP; ii++)
        {
            if ( IONS->at(ss).f3(ii) == 1 )
            {
                calculateRfTerms(ii,params,CS,fields,&IONS->at(ss));
            }

        } // OMP parallel for

    } // Species
}

void RF_Operator_TYP::calculatePowerPerUnitErf_AllSpecies(params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    double uEdot3 = 0;
    double DT = params->DT;

    for (int ss=0; ss<IONS->size();ss++)
    {
        double NCP = IONS->at(ss).NCP;
        int NSP    = IONS->at(ss).NSP;

        #pragma omp parallel default(none) shared(uEdot3, params, IONS, ss, CS, fields, std::cout) firstprivate(NSP,NCP,DT)
        {
            // Private variables:
            double uEdot3_private = 0;

            #pragma omp for
            for(int ii=0; ii<NSP; ii++)
            {
                if ( IONS->at(ss).f3(ii) == 1 )
                {
                    // Rf terms:
                    double udE3 = IONS->at(ss).udE3(ii);
                    double a_p  = IONS->at(ss).a_p(ii);

                    // Accumulate power:
                    uEdot3_private += (NCP/DT)*a_p*udE3;

                } // if

            } // omp for

            #pragma omp critical
            uEdot3 += uEdot3_private;

        } // omp parallel

    } // Species

    // Reduce over all MPI process
    MPI_AllreduceDouble(params,&uEdot3);

    // Assign to output:
    params->RF.uE3 = uEdot3;

    /*
    if (params->mpi.IS_PARTICLES_ROOT)
    {
        cout << "uEdot3 = " << params->RF.uE3*CS->energy/CS->time << endl;
    }
    */

}

void RF_Operator_TYP::ApplyRfOperator_AllSpecies( params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    // Calculate electric field:
    double E_rf = sqrt(params->RF.Prf/params->RF.uE3);

    /*
    if (params->mpi.IS_PARTICLES_ROOT)
    {
        cout << "E_RF = " << E_rf*CS->eField << endl;
    }
    */

    // Seed the random number generator:
    std::default_random_engine generator(params->mpi.MPI_DOMAIN_NUMBER+1);

    // Create uniform random number generator in [0,1]:
    std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

    for (int ss=0; ss<IONS->size();ss++)
    {
        int NSP = IONS->at(ss).NSP;
        int Ma  = IONS->at(ss).M;

        #pragma omp parallel default(none) shared(params, IONS, ss, CS, E_rf, NSP, Ma, cout, uniform_distribution) firstprivate(generator)
        {
            #pragma omp for
            for(int ii=0; ii<NSP; ii++)
            {
                if (IONS->at(ss).f3(ii) == 1)
                {
                    //  Particle states:
                    double vpar = IONS->at(ss).V_p(ii,0);
                    double vper = IONS->at(ss).V_p(ii,1);

                    // Sign of vpar:
                    double eps  = abs(vpar)/vpar;

                    // RF terms:
                    double mean_udKE_per = IONS->at(ss).udErf(ii);
                    double doppler       = IONS->at(ss).doppler(ii);

                    // Derived quantities:
                    double KE_par = Ma*vpar*vpar/2;
                    double KE_per = Ma*vper*vper/2;

                    // Calculate mean RF energy kick:
                    double mean_dKE_per = mean_udKE_per*pow(E_rf,2);

                    // Random number between 0 and 1:
                    double randomNumber = uniform_distribution(generator);
                    double Rm = (2*randomNumber - 1);

                    // Monte-Carlo operaton in kinetic energy:
                    double dKE_per = mean_dKE_per + Rm*sqrt(2*KE_per*mean_dKE_per);

                    // Total change in kinetic energy:
                    double dKE = dKE_per*(1 + doppler);

                    // Final perpendicular kinetic energy:
                    KE_per += dKE_per;

                    if (KE_per < 0)
                    {
                        cout << "KE_per is negative" << endl;
                    }

                    // Convert back to velocities:
                    vpar += (doppler/vpar)*(dKE_per/Ma);
                    vper = sqrt(2*KE_per/Ma);

                    // Output data:
                    IONS->at(ss).V_p(ii,0) = vpar;
                    IONS->at(ss).V_p(ii,1) = vper;

                    // Energy increments:
                    IONS->at(ss).dE3(ii) = dKE;

                    /*
                    if (params->mpi.IS_PARTICLES_ROOT)
                    {
                        cout << "mean_dKE_per [eV] = " << mean_dKE_per*CS->energy/F_E << endl;
                        cout << "dKE_per [eV] = " << dKE_per*CS->energy/F_E << endl;
                        cout << "dE3 [eV] = " << IONS->at(ss).dE3(ii)*CS->energy/F_E << endl;
                        cout << "dKE [eV] = " << dKE*CS->energy/F_E << endl;
                        cout << "vpar = " << vpar*CS->velocity << endl;

                        cout << "vper = " << vper*CS->velocity << endl;

                    }
                    */


                } // f3

            } // OMP for

        } // OMP parallel

    } // Species

}

void RF_Operator_TYP::calculateAbsorbedPower_AllSpecies(params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    double Edot3 = 0;
    double DT = params->DT;

    for (int ss=0; ss<IONS->size();ss++)
    {
        double NCP = IONS->at(ss).NCP;
        int NSP    = IONS->at(ss).NSP;

        #pragma omp parallel default(none) shared(Edot3, params, IONS, ss, CS, fields, std::cout) firstprivate(NSP,NCP,DT)
        {
            // Private variables:
            double Edot3_private = 0;

            #pragma omp for
            for(int ii=0; ii<NSP; ii++)
            {
                if ( IONS->at(ss).f3(ii) == 1 )
                {
                    // Rf terms:
                    double dE3 = IONS->at(ss).dE3(ii);
                    double a_p = IONS->at(ss).a_p(ii);

                    /*
                    if (params->mpi.IS_PARTICLES_ROOT)
                    {
                        cout << "a_p = " << a_p << endl;
                        cout << "dE3 [eV] = " << IONS->at(ss).dE3(ii)*CS->energy/F_E << endl;
                        cout << "vpar = " << IONS->at(ss).V_p(ii,0)*CS->velocity << endl;
                        cout << "vper = " << IONS->at(ss).V_p(ii,1)*CS->velocity << endl;
                    }
                    */

                    // Accumulate power:
                    Edot3_private += (NCP/DT)*a_p*dE3;

                    // Clear flags:
                    IONS->at(ss).f3(ii)  = 0;
                    IONS->at(ss).dE3(ii) = 0;

                } // if

            } // omp for

            #pragma omp critical
            {
                Edot3 += Edot3_private;
            }

        } // omp parallel

    } // Species

    // Reduce over all MPI process
    MPI_AllreduceDouble(params,&Edot3);

    // Assign to output:
    params->RF.E3 = Edot3;


    if (params->mpi.IS_PARTICLES_ROOT)
    {
        //cout << (params->RF.E3*CS->energy/CS->time)/1000 << endl;
    }

}

void RF_Operator_TYP::ApplyRfHeating_AllSpecies(params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        // Calculate Resonance number:
        calculateResNum_AllSpecies(params,CS,fields,IONS);

        // Check resonance condition and flag:
        checkResNumAndFlag_AllSpecies(params,CS,fields,IONS);

        // calculate RF terms and unit kick:
        calculateRfTerms_AllSpecies(params,CS,fields,IONS);

        // Calculate RF power per unit electric field over all species:
        calculatePowerPerUnitErf_AllSpecies(params,CS,fields,IONS);

        // Apply RF heating to all allSpecies:
        ApplyRfOperator_AllSpecies(params,CS,fields,IONS);

        // Calculate the absorbed RF power over all species:
        calculateAbsorbedPower_AllSpecies(params,CS,fields,IONS);

    } // PARTICLE MPI
}
