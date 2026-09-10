#include "rfOperator.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>

#ifndef HAS_STD_BESSEL
#include <boost/math/special_functions/bessel.hpp>
#endif

namespace
{
RF_SPECIES_TYP& rfConfigForSpecies(params_TYP * params, const ionSpecies_TYP & species)
{
    if (species.Z < 0.0)
    {
        return params->RF.electrons;
    }
    return params->RF.ions;
}

const RF_SPECIES_TYP& rfConfigForSpecies(const params_TYP * params, const ionSpecies_TYP & species)
{
    if (species.Z < 0.0)
    {
        return params->RF.electrons;
    }
    return params->RF.ions;
}

bool rfHeatingEnabledForSpecies(const params_TYP * params, const ionSpecies_TYP & species)
{
    if (species.Z > 0.0)
    {
        return params->RF.heatIons == 1;
    }
    if (species.Z < 0.0)
    {
        return params->RF.heatElectrons == 1;
    }
    return false;
}

bool rfHeatingActiveForSpecies(const params_TYP * params, const CS_TYP * CS, const ionSpecies_TYP & species)
{
    if (!rfHeatingEnabledForSpecies(params, species))
    {
        return false;
    }
    const RF_SPECIES_TYP& rf = rfConfigForSpecies(params, species);
    return params->currentTime >= rf.t_ON*CS->time && params->currentTime <= rf.t_OFF*CS->time;
}

bool relativisticElectronsEnabledForSpecies(const params_TYP * params, const ionSpecies_TYP & species)
{
    return params->SW.relativisticElectrons == 1 && species.Z < 0.0;
}

double perpendicularSpeedForRf(const params_TYP * params, const ionSpecies_TYP & species, int ii)
{
    if (params->advanceParticleMethod == PARTICLE_PUSH_BORIS_FULL_ORBIT && species.V_p.n_cols > 2)
    {
        return hypot(species.V_p(ii,1), species.V_p(ii,2));
    }
    return fabs(species.V_p(ii,1));
}

double safeBesselJ(int order, double argument)
{
    if (order < 0 || !std::isfinite(argument) || argument < 0.0)
    {
        return 0.0;
    }

    try
    {
        return CYL_BESSEL_J(order, argument);
    }
    catch (const std::exception&)
    {
        return 0.0;
    }
}

double gammaFromSpeed(double speed)
{
    const double c = std::max(F_C_DS, double_zero);
    double beta2 = speed*speed/(c*c);
    beta2 = std::max(0.0, std::min(beta2, 1.0 - 1.0e-12));
    return 1.0/sqrt(1.0 - beta2);
}

double gammaFromVelocity(double vpar, double vper)
{
    return gammaFromSpeed(hypot(vpar, vper));
}

double kineticEnergyFromSpeed(double mass, double speed, bool relativistic)
{
    if (!relativistic)
    {
        return 0.5*mass*speed*speed;
    }
    const double gamma = gammaFromSpeed(speed);
    return (gamma - 1.0)*mass*F_C_DS*F_C_DS;
}

double speedFromKineticEnergy(double mass, double kineticEnergy, bool relativistic)
{
    if (kineticEnergy <= 0.0)
    {
        return 0.0;
    }
    if (!relativistic)
    {
        return sqrt(2.0*kineticEnergy/mass);
    }

    const double restEnergy = mass*F_C_DS*F_C_DS;
    const double gamma = 1.0 + kineticEnergy/restEnergy;
    const double beta2 = std::max(0.0, std::min(1.0 - 1.0e-12, 1.0 - 1.0/(gamma*gamma)));
    return F_C_DS*sqrt(beta2);
}

double perpendicularKineticEnergy(double mass, double vpar, double vper, bool relativistic)
{
    if (!relativistic)
    {
        return 0.5*mass*vper*vper;
    }

    const double speed2 = vpar*vpar + vper*vper;
    if (speed2 <= double_zero)
    {
        return 0.0;
    }
    const double totalEnergy = kineticEnergyFromSpeed(mass, sqrt(speed2), true);
    return totalEnergy*(vper*vper/speed2);
}
}

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
    double vper = perpendicularSpeedForRf(params, *IONS, ii);
    double Bp   = fabs(IONS->BX_p(ii));
    if (!std::isfinite(Bp) || Bp <= double_zero || !std::isfinite(vpar) || !std::isfinite(vper))
    {
        IONS->resNum(ii) = std::numeric_limits<double>::infinity();
        return;
    }
    double gamma = relativisticElectronsEnabledForSpecies(params, *IONS) ? gammaFromVelocity(vpar, vper) : 1.0;
    if (!std::isfinite(gamma) || gamma <= double_zero)
    {
        IONS->resNum(ii) = std::numeric_limits<double>::infinity();
        return;
    }
    double wcp   = abs(Q)*Bp/(gamma*Ma);

    // RF paramters:
    const RF_SPECIES_TYP& rf = rfConfigForSpecies(params, *IONS);
    int n       = rf.n_harmonic;
    double kpar = rf.kpar;
    double f    = rf.freq;
    double wrf  = 2*M_PI*f;

    // Cyclotron resonance number:
    IONS->resNum(ii) = wrf - kpar*vpar - n*wcp;
}

void RF_Operator_TYP::calculateResNum_AllSpecies(params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
        for (int ss=0; ss<IONS->size();ss++)
        {
            if (!rfHeatingActiveForSpecies(params, CS, IONS->at(ss)))
            {
                IONS->at(ss).f3.zeros();
                IONS->at(ss).dE3.zeros();
                IONS->at(ss).udErf.zeros();
                IONS->at(ss).udE3.zeros();
                continue;
            }

            int NSP   = IONS->at(ss).NSP;

            #pragma omp parallel for default(none) shared(params, IONS, ss, CS, fields, std::cout) firstprivate(NSP)
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
        if (!rfHeatingActiveForSpecies(params, CS, IONS->at(ss)))
        {
            continue;
        }

        int NSP   = IONS->at(ss).NSP;

        #pragma omp parallel for default(none) shared(params, IONS, ss, CS, fields, std::cout) firstprivate(NSP)
        for(int ii=0; ii<NSP; ii++)
        {
            // Parameters needed to check resonance condition:
            double xp = IONS->at(ss).X_p(ii);
            const RF_SPECIES_TYP& rf = rfConfigForSpecies(params, IONS->at(ss));
            double x1 = rf.x1;
            double x2 = rf.x2;
            double resNum  = IONS->at(ss).resNum(ii);
            double resNum_ = IONS->at(ss).resNum_(ii);

            bool isResonant = false;
            if (rf.resonanceMode == RF_RESONANCE_FORTRAN_WINDOW)
            {
                isResonant = resNum < 0.0;
            }
            else
            {
                // true: argument negative; false: otherwise (positive or zero)
                isResonant = signbit(resNum*resNum_);
            }

            // Flag particles in resonance:
            if (isResonant && (xp > x1) && (xp < x2))
            {
                IONS->at(ss).f3(ii) = 1;
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
    double Z  = IONS->Z;

    // RF paramters:
    const RF_SPECIES_TYP& rf = rfConfigForSpecies(params, *IONS);
    int n         = rf.n_harmonic;
    double kper   = rf.kper;
    double kpar   = rf.kpar;
    double tau_rf = 0;
    double rL     = 0;
    double flr    = 0;
    double J_nm1  = 0;
    double J_np1  = 0;
    double mean_dKE_per = 0;

    // Particle states:
    double vpar = IONS->V_p(ii,0);
    double vper = perpendicularSpeedForRf(params, *IONS, ii);
    const bool relativistic = relativisticElectronsEnabledForSpecies(params, *IONS);
    const double gamma = relativistic ? gammaFromVelocity(vpar, vper) : 1.0;

    // Particle-defined fields:
    double Bp   = fabs(IONS->BX_p(ii));
    double dBp  = IONS->dBX_p(ii);
    double ddBp = IONS->ddBX_p(ii);
    double Ep   = IONS->EX_p(ii);
    if (!std::isfinite(Bp) || Bp <= double_zero ||
        !std::isfinite(dBp) || !std::isfinite(ddBp) || !std::isfinite(Ep) ||
        !std::isfinite(vpar) || !std::isfinite(vper) ||
        !std::isfinite(gamma) || gamma <= double_zero || Ma <= double_zero ||
        n <= 0)
    {
        IONS->udErf(ii) = 0.0;
        IONS->doppler(ii) = 0.0;
        IONS->udE3(ii) = 0.0;
        return;
    }

    // Derived quantities:
    double Omega   = abs(Q)*Bp/(gamma*Ma);
    double dOmega  = abs(Q)*dBp/(gamma*Ma);
    double ddOmega = abs(Q)*ddBp/(gamma*Ma);
    double qOverMass = Q/(gamma*Ma);
    if (!std::isfinite(Omega) || Omega <= double_zero ||
        !std::isfinite(dOmega) || !std::isfinite(ddOmega) || !std::isfinite(qOverMass))
    {
        IONS->udErf(ii) = 0.0;
        IONS->doppler(ii) = 0.0;
        IONS->udE3(ii) = 0.0;
        return;
    }

    // Calculate the first and second time derivative of Omega:
    double Omega_dot  = vpar*dOmega;
    double Omega_ddot = pow(vpar,2)*ddOmega  - pow(vper,2)*pow(dOmega,2)/(2.*Omega)  +  qOverMass*Ep*dOmega;
    if (!std::isfinite(Omega_dot) || !std::isfinite(Omega_ddot))
    {
        IONS->udErf(ii) = 0.0;
        IONS->doppler(ii) = 0.0;
        IONS->udE3(ii) = 0.0;
        return;
    }

    // Calculate the interaction time:
    if ( pow(n*Omega_ddot,2) > 4.8175*abs(pow(n*Omega_dot,3)) )
    {
        // tau_b
        // Approximate Ai(x) ~ 0.3833
        if (abs(n*Omega_ddot) <= double_zero)
        {
            IONS->udErf(ii) = 0.0;
            IONS->doppler(ii) = 0.0;
            IONS->udE3(ii) = 0.0;
            return;
        }
        tau_rf = (2*M_PI)*pow(abs(2/(n*Omega_ddot)),1.0/3.0)*0.3833;
    }
    else
    {
        // tau_a
        if (abs(n*Omega_dot) <= double_zero)
        {
            IONS->udErf(ii) = 0.0;
            IONS->doppler(ii) = 0.0;
            IONS->udE3(ii) = 0.0;
            return;
        }
        tau_rf = sqrt(2*M_PI/abs(n*Omega_dot));
    }
    if (!std::isfinite(tau_rf) || tau_rf <= double_zero)
    {
        IONS->udErf(ii) = 0.0;
        IONS->doppler(ii) = 0.0;
        IONS->udE3(ii) = 0.0;
        return;
    }

    // Calculate bessel terms:
    rL  = vper/Omega;
    flr = fabs(kper)*rL;
    if (!std::isfinite(flr) || flr < 0.0)
    {
        IONS->udErf(ii) = 0.0;
        IONS->doppler(ii) = 0.0;
        IONS->udE3(ii) = 0.0;
        return;
    }
    J_nm1 = safeBesselJ(n - 1, flr);
    J_np1 = safeBesselJ(n + 1, flr);

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
        mean_dKE_per = 0.5*(pow(e,2)/(gamma*Ma))*pow(abs(E_p*J_nm1 + E_m*J_np1)*tau_rf,2); // [J] normalized energy
    }
    if (Z < 0) // Negative particles
    {
        double E_m = 1;
        double E_p = 0;
        mean_dKE_per = 0.5*(pow(e,2)/(gamma*Ma))*pow(abs(E_m*J_nm1 + E_p*J_np1)*tau_rf,2); // [J] normalized energy
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
        if (!rfHeatingActiveForSpecies(params, CS, IONS->at(ss)))
        {
            continue;
        }

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
    double ionUEdot3 = 0;
    double electronUEdot3 = 0;
    double DT = params->DT;

    for (int ss=0; ss<IONS->size();ss++)
    {
        if (!rfHeatingActiveForSpecies(params, CS, IONS->at(ss)))
        {
            continue;
        }

        double& speciesUEdot3 = (IONS->at(ss).Z < 0.0) ? electronUEdot3 : ionUEdot3;
        double NCP = IONS->at(ss).NCP;
        int NSP    = IONS->at(ss).NSP;

        #pragma omp parallel default(none) shared(speciesUEdot3, params, IONS, ss, CS, fields, std::cout) firstprivate(NSP,NCP,DT)
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
            speciesUEdot3 += uEdot3_private;

        } // omp parallel

    } // Species

    // Reduce over all MPI process
    MPI_AllreduceDouble(params,&ionUEdot3);
    MPI_AllreduceDouble(params,&electronUEdot3);

    // Assign to output:
    params->RF.ions.uE3 = ionUEdot3;
    params->RF.electrons.uE3 = electronUEdot3;
    params->RF.uE3 = ionUEdot3 + electronUEdot3;

}

void RF_Operator_TYP::ApplyRfOperator_AllSpecies( params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    // Seed the random number generator:
    std::default_random_engine generator(params->mpi.MPI_DOMAIN_NUMBER+1);

    // Create uniform random number generator in [0,1]:
    std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

    for (int ss=0; ss<IONS->size();ss++)
    {
        if (!rfHeatingActiveForSpecies(params, CS, IONS->at(ss)))
        {
            continue;
        }

        RF_SPECIES_TYP& rf = rfConfigForSpecies(params, IONS->at(ss));
        double E_rf = rf.eFieldAmplitude;
        if (rf.eFieldMode == RF_EFIELD_POWER_BALANCE)
        {
            if ((rf.uE3 <= double_zero) || !std::isfinite(rf.uE3))
            {
                rf.Erf = 0.0;
                params->RF.Erf = 0.0;
                continue;
            }
            E_rf = sqrt(std::max(0.0, rf.Prf/rf.uE3));
        }
        rf.Erf = E_rf;
        params->RF.Erf = E_rf;

        const bool relativistic = relativisticElectronsEnabledForSpecies(params, IONS->at(ss));
        int NSP = IONS->at(ss).NSP;
        double Ma  = IONS->at(ss).M;
        double maxEnergyGainFraction = rf.maxEnergyGainFraction;
        double maxParticleEnergy = rf.maxParticleEnergy;
        double maxVelocityFractionC = rf.maxVelocityFractionC;

        #pragma omp parallel default(none) shared(params, IONS, ss, CS, E_rf, NSP, Ma, cout, uniform_distribution, maxEnergyGainFraction, maxParticleEnergy, maxVelocityFractionC, F_C_DS) firstprivate(generator, relativistic)
        {
            #pragma omp for
            for(int ii=0; ii<NSP; ii++)
            {
                if (IONS->at(ss).f3(ii) == 1)
                {
                    //  Particle states:
                    double vpar = IONS->at(ss).V_p(ii,0);
                    double vper = perpendicularSpeedForRf(params, IONS->at(ss), ii);

                    // Sign of vpar:
                    double eps  = (vpar >= 0.0) ? 1.0 : -1.0;

                    // RF terms:
                    double mean_udKE_per = IONS->at(ss).udErf(ii);
                    double doppler       = IONS->at(ss).doppler(ii);

                    // Derived quantities:
                    double KE_per = perpendicularKineticEnergy(Ma, vpar, vper, relativistic);
                    const double KE_total_before = kineticEnergyFromSpeed(Ma, hypot(vpar, vper), relativistic);

                    // Calculate mean RF energy kick:
                    double mean_dKE_per = mean_udKE_per*pow(E_rf,2);

                    // Random number between 0 and 1:
                    double randomNumber = uniform_distribution(generator);
                    double Rm = (2*randomNumber - 1);

                    // Monte-Carlo operaton in kinetic energy:
                    double dKE_per = mean_dKE_per + Rm*sqrt(2*KE_per*mean_dKE_per);
                    if (maxEnergyGainFraction > 0.0)
                    {
                        const double referenceEnergy = std::max(KE_per, params->f_IC.Te);
                        const double maxPositiveKick = maxEnergyGainFraction*referenceEnergy;
                        dKE_per = std::min(dKE_per, maxPositiveKick);
                    }
                    dKE_per = std::max(dKE_per, -0.95*KE_per);

                    // Total change in kinetic energy:
                    double dKE = dKE_per*(1 + doppler);

                    if (relativistic)
                    {
                        double targetTotalKE = std::max(0.0, KE_total_before + dKE);
                        double targetPerpKE = std::max(0.0, KE_per + dKE_per);
                        if (targetTotalKE < targetPerpKE)
                        {
                            targetTotalKE = targetPerpKE;
                        }

                        if (maxParticleEnergy > 0.0)
                        {
                            targetTotalKE = std::min(targetTotalKE, maxParticleEnergy);
                        }
                        if (maxVelocityFractionC > 0.0)
                        {
                            const double maxSpeed = std::min(maxVelocityFractionC, F_C_DS*(1.0 - 1.0e-9));
                            const double maxSpeedEnergy = kineticEnergyFromSpeed(Ma, maxSpeed, true);
                            targetTotalKE = std::min(targetTotalKE, maxSpeedEnergy);
                        }

                        targetPerpKE = std::min(targetPerpKE, targetTotalKE);
                        const double newSpeed = speedFromKineticEnergy(Ma, targetTotalKE, true);
                        const double perpFraction = (targetTotalKE > double_zero) ? std::max(0.0, std::min(1.0, targetPerpKE/targetTotalKE)) : 0.0;
                        vper = newSpeed*sqrt(perpFraction);
                        vpar = eps*newSpeed*sqrt(std::max(0.0, 1.0 - perpFraction));
                        dKE = targetTotalKE - KE_total_before;
                    }
                    else
                    {
                        const double KE_par = std::max(0.0, KE_total_before - KE_per);
                        double targetPerpKE = KE_per + dKE_per;
                        double targetParKE = KE_par + doppler*dKE_per;

                        if (targetPerpKE < 0)
                        {
                            cout << "KE_per is negative" << endl;
                            targetPerpKE = 0.0;
                        }
                        if (targetParKE < 0)
                        {
                            targetParKE = 0.0;
                        }

                        double totalKE = targetPerpKE + targetParKE;
                        double cappedKE = totalKE;
                        if (maxParticleEnergy > 0.0)
                        {
                            cappedKE = std::min(cappedKE, maxParticleEnergy);
                        }
                        if (maxVelocityFractionC > 0.0)
                        {
                            const double maxSpeedEnergy = 0.5*Ma*maxVelocityFractionC*maxVelocityFractionC;
                            cappedKE = std::min(cappedKE, maxSpeedEnergy);
                        }
                        if (cappedKE < totalKE && totalKE > double_zero)
                        {
                            const double scaleEnergy = cappedKE/totalKE;
                            targetPerpKE *= scaleEnergy;
                            targetParKE *= scaleEnergy;
                            totalKE = cappedKE;
                            dKE = cappedKE - KE_total_before;
                        }

                        vpar = eps*sqrt(2.0*targetParKE/Ma);
                        vper = sqrt(2.0*targetPerpKE/Ma);
                    }

                    // Output data:
                    IONS->at(ss).V_p(ii,0) = vpar;
                    if (params->advanceParticleMethod == PARTICLE_PUSH_BORIS_FULL_ORBIT && IONS->at(ss).V_p.n_cols > 2)
                    {
                        const double oldVper = hypot(IONS->at(ss).V_p(ii,1), IONS->at(ss).V_p(ii,2));
                        if (oldVper > double_zero)
                        {
                            const double scale = vper/oldVper;
                            IONS->at(ss).V_p(ii,1) *= scale;
                            IONS->at(ss).V_p(ii,2) *= scale;
                        }
                        else
                        {
                            IONS->at(ss).V_p(ii,1) = vper;
                            IONS->at(ss).V_p(ii,2) = 0.0;
                        }
                    }
                    else
                    {
                        IONS->at(ss).V_p(ii,1) = vper;
                    }

                    // Energy increments:
                    IONS->at(ss).dE3(ii) = dKE;

                } // f3

            } // OMP for

        } // OMP parallel

    } // Species

}

void RF_Operator_TYP::calculateAbsorbedPower_AllSpecies(params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    double ionEdot3 = 0;
    double electronEdot3 = 0;
    double DT = params->DT;

    for (int ss=0; ss<IONS->size();ss++)
    {
        if (!rfHeatingActiveForSpecies(params, CS, IONS->at(ss)))
        {
            continue;
        }

        double& speciesEdot3 = (IONS->at(ss).Z < 0.0) ? electronEdot3 : ionEdot3;
        double NCP = IONS->at(ss).NCP;
        int NSP    = IONS->at(ss).NSP;

        #pragma omp parallel default(none) shared(speciesEdot3, params, IONS, ss, CS, fields, std::cout) firstprivate(NSP,NCP,DT)
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

                    // Accumulate power:
                    Edot3_private += (NCP/DT)*a_p*dE3;

                    // Clear flags:
                    IONS->at(ss).f3(ii)  = 0;
                    IONS->at(ss).dE3(ii) = 0;

                } // if

            } // omp for

            #pragma omp critical
            {
                speciesEdot3 += Edot3_private;
            }

        } // omp parallel

    } // Species

    // Reduce over all MPI process
    MPI_AllreduceDouble(params,&ionEdot3);
    MPI_AllreduceDouble(params,&electronEdot3);

    // Assign to output:
    params->RF.ions.E3 = ionEdot3;
    params->RF.electrons.E3 = electronEdot3;
    params->RF.E3 = ionEdot3 + electronEdot3;

    /*
    if (params->mpi.IS_PARTICLES_ROOT)
    {
        cout << (params->RF.E3*CS->energy/CS->time)/1000 << endl;
    }
    */

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

        for (ionSpecies_TYP &ion : *IONS)
        {
            if (!rfHeatingActiveForSpecies(params, CS, ion))
            {
                continue;
            }
            const RF_SPECIES_TYP& rf = rfConfigForSpecies(params, ion);
            if (rf.eFieldMode == RF_EFIELD_POWER_BALANCE &&
                ((rf.uE3 <= double_zero) || !std::isfinite(rf.uE3)))
            {
                ion.f3.zeros();
                ion.dE3.zeros();
            }
        }

        // Apply RF heating to all allSpecies:
        ApplyRfOperator_AllSpecies(params,CS,fields,IONS);

        // Calculate the absorbed RF power over all species:
        calculateAbsorbedPower_AllSpecies(params,CS,fields,IONS);

    } // PARTICLE MPI
}
