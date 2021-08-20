#include "rfOperator.h"

RF_Operator_TYP::RF_Operator_TYP(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
}

void RF_Operator_TYP::cyclotronResonanceNumber(int ii, const params_TYP * params, CS_TYP * CS, fields_TYP * fields, ionSpecies_TYP * IONS)
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
    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
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
                cyclotronResonanceNumber(ii,params,CS,fields,&IONS->at(ss));

            } // particles
        } // Species
    } // PARTICLE MPI
}



void RF_Operator_TYP::ApplyRfHeating(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    // Calculate Resonance number:
    calculateResNum_AllSpecies(params,CS,fields,IONS);
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
