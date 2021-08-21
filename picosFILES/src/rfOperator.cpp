#include "rfOperator.h"

RF_Operator_TYP::RF_Operator_TYP(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {

    } // PARTICLE MPI
}

void RF_Operator_TYP::calculateResNum(int ii, const params_TYP * params, CS_TYP * CS, fields_TYP * fields, ionSpecies_TYP * IONS)
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

void RF_Operator_TYP::calculateResNum_AllSpecies(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
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

void RF_Operator_TYP::checkResNumAndFlag_AllSpecies(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
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

void RF_Operator_TYP::calculateRfTerms(int ii, const params_TYP * params, CS_TYP * CS, fields_TYP * fields, ionSpecies_TYP * IONS)
{
    // Ion parameters:
    double Ma = IONS->M;
    double Q  = IONS->Q;
    double e  = F_E_DS;

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
    double KE_par   = 0.5*Ma*vper*vper;
    double KE_per   = 0.5*Ma*vper*vper;
    double Omega    = abs(Q)*Bp/Ma;
    double nOmega   = n*Omega;
    double ndOmega  = n*abs(Q)*dBp/Ma;
    double nddOmega = n*abs(Q)*ddBp/Ma;

    // Calculate the first and second time derivative of Omega:
    double nOmega_dot  = vpar*ndOmega;
    double nOmega_ddot = pow(vpar,2)*nddOmega  - pow(vper,2)*pow(ndOmega,2)/(2.*nOmega) + Q*Ep*ndOmega/Ma;

    // Calculate the interaction time:
    if ( pow(nOmega_ddot,2) > 4.8175*abs(pow(nOmega_dot,3)) )
    {
        // tau_b
        // Approximate Ai(x) ~ 0.3833
        //tau_rf = (2*M_PI)*(abs(2/nOmega_ddot)**(1/3.))*0.3833;
        tau_rf = (2*M_PI)*pow(abs(2/nOmega_ddot),1/3)*0.3833;

        // *******************
        // Need to check that the above fractional power is correct!!!
        // ***************


        cout << "tau_b = " << tau_rf*CS->time << "[s]" << endl;

    }
    else
    {
        // tau_a
        tau_rf = sqrt(2*M_PI/abs(nOmega_dot));
        cout << "tau_a = " << tau_rf*CS->time << "[s]" << endl;

    }

    // Calculate bessel terms:
    rL  = vper/Omega;
    flr = kper*rL;
    J_nm1 = std::cyl_bessel_j(n - 1,flr);
    J_np1 = std::cyl_bessel_j(n + 1,flr);

    // Calculate the mean RF kick per unit electric field squared:
    // ********************
    // This applies to ions only!!
    // Need to generalize to accept electrons
    // *********************
    double E_m = 0;
    double E_p = 1;
    mean_dKE_per = 0.5*(e/Ma)*pow(abs(E_p*J_nm1 + E_m*J_np1)*tau_rf,2);

    // Populate output:
    //IONS->udErf(ii)   = mean_dKE_per;
    //IONS->doppler(ii) = kpar*vpar/nOmega;
    //IONS->udE3(ii)    =IONS->udErf(ii)*(1. + IONS->doppler(ii));

    // *********************
    // Need to create the above variables in IONS
    // **********************

}

void RF_Operator_TYP::calculateRfTerms_AllSpecies(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
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

        } // particles
    } // Species
}

void RF_Operator_TYP::ApplyRfHeating(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        // Calculate Resonance number:
        calculateResNum_AllSpecies(params,CS,fields,IONS);

        // Check resonance condition and flag:
        checkResNumAndFlag_AllSpecies(params,CS,fields,IONS);

        // calculate RF terms and unit kick:
        calculateRfTerms_AllSpecies(params,CS,fields,IONS);

    } // PARTICLE MPI
}



// Calculate resNum in RF operator constructor
// rf_operator.ApplyRfHeating()
// {
//    calciulateResNum_AllSpecies()
//   For all species and for all particles (OMP)
//   {
//      calculateResNumber()
//      {
//           IONS->at(ss).resNum_(ii) = IONS->at(ss).resNum(ii);
//          calculateResNum(params,fields,IONS);
//      }
//    }


//   calculateRFterms_AllSpecies()
//   For all species and for all particles (OMP)
//   {
//      checkResonanceCondition(params,fields,IONS)
//      {
//          Compare resNum_ and resNum, have they changed sign?
//          {
//              IONS->at(ss).f3(ii) = 1;
//              calculateRFterms();
//              {
//                  Apply unit RF kick
//                  Calculate interatction time, besselterms, dopplerterms
//                  IONS->at(ss).u_mean_dErf(ii);
//                  IONS->at(ss).udE3(ii);
//              }
//           }
//      }
//
//   calculatePowerPerUnitErf();
//   {
//   global_uEdot3 - 0;
//   For all species
//   {
//      local_uEdot3 = 0;
//      For all particles (OMP)
//      {
//          local_uEdot3 += udE3(ii)*otherTerms like alpha/dt
//      }
//      MPI_AllreduceDouble(local_uEdot3)
//   }
//   }
//
//  ApplyRF_AllSpecies(params,fields,IONS,global_uEdot3)
// {
//  Calculate E_rf based on global_uEdot3 and params->P_RF
//  loop over all species
//      loop over all particles(OMP)
//      {
//          if (f3(ii) == 1)
//          {
//              Apply monte carlo RF operatio bsded om E_rf and produce Eperp and Eparp
//              Comvert to Vpar and Vper
//              IONS->at(ss).dE3 = dKE;
//          }
//      }
// }
//
// calculateAbsorbedPower_AllSpecies(params,fields,IONS)
// {
//      Edot_3_global = 0;
//      Loop over all Species
//      {
//         loop over all PARTICLES
//          {
//              Edot_3_local += dE3*terms
//          }
//          MPI_AllreduceDouble(Edot_3_global,Edot_3_local)
//      }
//      params->RF.Edot3 =  Edot_3_global;
//  }
