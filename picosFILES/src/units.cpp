#include "units.h"

void units_TYP::defineCharacteristicScalesAndBcast(params_TYP * params, vector<ionSpecies_TYP> * IONS, CS_TYP * CS)
{
    MPI_Barrier(MPI_COMM_WORLD);

    if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        defineCharacteristicScales(params, IONS, CS);
    }

    broadcastCharacteristicScales(params, CS);
}


void units_TYP::defineCharacteristicScales(params_TYP * params, vector<ionSpecies_TYP> * IONS, CS_TYP * CS)
{
	// The definition of the characteristic quantities is based on:
	// D Winske and N Omidi, Hybrid codes.
	// All the quantities below have units (SI).

	cout << endl << "* * * * * * * * * * * * DEFINING CHARACTERISTIC SCALES IN SIMULATION * * * * * * * * * * * * * * * * * *" << endl;

	for (int ii=0; ii<params->numberOfParticleSpecies; ii++)
    {
	    CS->mass += IONS->at(ii).M;
	    CS->charge += fabs(IONS->at(ii).Q);
	}

	CS->mass   /= params->numberOfParticleSpecies;
	CS->charge /= params->numberOfParticleSpecies;
	CS->density = params->CV.ne;

	double characteristicPlasmaFrequency(0);//Background ion-plasma frequency.
	characteristicPlasmaFrequency = sqrt( CS->density*CS->charge*CS->charge/(CS->mass*F_EPSILON) );

	CS->time = 1/characteristicPlasmaFrequency;
	CS->velocity = F_C;
	CS->momentum = CS->mass*CS->velocity;
	CS->length = CS->velocity*CS->time;
	CS->volume = CS->length*CS->length*CS->length;
	CS->eField = ( CS->mass*CS->velocity )/( CS->charge*CS->time );
	CS->bField = CS->eField/CS->velocity; // CS->mass/( CS->charge*CS->time );
	CS->temperature = CS->mass*CS->velocity*CS->velocity/F_KB;
	CS->pressure = CS->bField*CS->velocity*CS->velocity*CS->charge*CS->density*CS->time;
	CS->resistivity = CS->bField/(CS->charge*CS->density);
	CS->magneticMoment = CS->mass*CS->velocity*CS->velocity/CS->bField;
	CS->vacuumPermittivity = (pow(CS->length*CS->charge,2)*CS->density)/(CS->mass*pow(CS->velocity,2));
	CS->vacuumPermeability = CS->mass/( CS->density*pow(CS->charge*CS->velocity*CS->time,2) );
    CS->energy = CS->mass*pow(CS->length/CS->time,2);

	if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
		cout << "+ Average mass: " << scientific << CS->mass << fixed << " kg" << endl;
		cout << "+ Average charge: " << scientific << CS->charge << fixed << " C" << endl;
		cout << "+ Density: " << scientific << CS->density << fixed << " m^(-3)" << endl;
		cout << "+ Time: " << scientific << CS->time << fixed << " s" << endl;
		cout << "+ Plasma frequency: " << scientific << characteristicPlasmaFrequency << fixed << " s" << endl;
		cout << "+ Velocity: " << scientific << CS->velocity << fixed << " m/s" << endl;
		cout << "+ Length: " << scientific << CS->length << fixed << " m" << endl;
		cout << "+ Electric field intensity: " << scientific << CS->eField << fixed << " V/m" << endl;
		cout << "+ Magnetic field intensity: " << scientific << CS->bField << fixed << " T" << endl;
		cout << "+ Pressure: " << scientific << CS->pressure << fixed << " Pa" << endl;
		cout << "+ Temperature: " << scientific << CS->temperature << fixed << " K" << endl;
		cout << "+ Magnetic moment: " << scientific << CS->magneticMoment << fixed << " A*m^2" << endl;
		cout << "+ Resistivity: " << scientific << CS->magneticMoment << fixed << " Ohms*m" << endl;
		cout << "+ Vacuum permittivity: " << scientific << CS->vacuumPermittivity << fixed << "" << endl;
		cout << "+ Vacuum permeability: " << scientific << CS->vacuumPermeability << fixed << "" << endl;
		cout << endl << "* * * * * * * * * * * * CHARACTERISTIC SCALES IN SIMULATION DEFINED  * * * * * * * * * * * * * * * * * *" << endl;
	}
}



void units_TYP::broadcastCharacteristicScales(params_TYP * params, CS_TYP * CS)
{
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&CS->time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->velocity, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->momentum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->length, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Bcast(&CS->volume, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->mass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->charge, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->density, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->eField, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->bField, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->pressure, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->temperature, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->magneticMoment, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->resistivity, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->vacuumPermeability, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->vacuumPermittivity, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&CS->energy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}



void units_TYP::calculateFundamentalScalesAndBcast(params_TYP * params, vector<ionSpecies_TYP> * IONS, FS_TYP * FS)
{
    MPI_Barrier(MPI_COMM_WORLD);

    if(params->mpi.MPI_DOMAIN_NUMBER == 0)
	{
        calculateFundamentalScales(params, IONS, FS);
    }

    broadcastFundamentalScales(params, FS);

    // It is assumed that species 0 is the majority species
    params->ionLarmorRadius = FS->ionGyroRadius[0];
    params->ionSkinDepth = FS->ionSkinDepth[0];
    params->ionGyroPeriod = FS->ionGyroPeriod[0];
}



void units_TYP::calculateFundamentalScales(params_TYP * params, vector<ionSpecies_TYP> * IONS, FS_TYP * FS)
{

	cout << endl << "* * * * * * * * * * * * CALCULATING FUNDAMENTAL SCALES IN SIMULATION * * * * * * * * * * * * * * * * * *" << endl;

	FS->electronSkinDepth  = F_C/sqrt( params->CV.ne*F_E*F_E/(F_EPSILON*F_ME) );
	FS->electronGyroPeriod = 2.0*M_PI/(F_E*params->CV.B/F_ME);
	FS->electronGyroRadius = sqrt(2.0*F_KB*params->CV.Te/F_ME)/(F_E*params->CV.B/F_ME);

	cout << " + Electron gyro-period: " << scientific << FS->electronGyroPeriod << fixed << " s" << endl;
	cout << " + Electron skin depth: " << scientific << FS->electronSkinDepth << fixed << " m" << endl;
	cout << " + Electron gyro-radius: " << scientific << FS->electronGyroRadius << fixed << " m" << endl;
	cout << endl;

	for(int ss=0; ss<params->numberOfParticleSpecies; ss++){
		FS->ionGyroPeriod[ss] = 2.0*M_PI/IONS->at(ss).Wc;
		FS->ionSkinDepth[ss]  = F_C/IONS->at(ss).Wp;
		FS->ionGyroRadius[ss] = IONS->at(ss).LarmorRadius;

		cout << "ION SPECIES: " << ss << endl;
		cout << " + Ion gyro-period: " << scientific << FS->ionGyroPeriod[ss] << fixed << " s" << endl;
		cout << " + Ion skin depth: " << scientific << FS->ionSkinDepth[ss] << fixed << " m" << endl;
		cout << " + Ion gyro-radius: " << scientific << FS->ionGyroRadius[ss] << fixed << " m" << endl;
		cout << endl;
	}

	cout << "* * * * * * * * * * * * FUNDAMENTAL SCALES IN SIMULATION CALCULATED  * * * * * * * * * * * * * * * * * *" << endl;
}



void units_TYP::broadcastFundamentalScales(params_TYP * params, FS_TYP * FS)
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Electron skin depth
    MPI_Bcast(&FS->electronSkinDepth, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Electron gyro-period
    MPI_Bcast(&FS->electronGyroPeriod, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Electron gyro-radius
    MPI_Bcast(&FS->electronGyroRadius, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Ion skin depth
    MPI_Bcast(FS->ionSkinDepth, params->numberOfParticleSpecies, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Ion gyro-period
    MPI_Bcast(FS->ionGyroPeriod, params->numberOfParticleSpecies, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Ion gyro-radius
    MPI_Bcast(FS->ionGyroRadius, params->numberOfParticleSpecies, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}


void units_TYP::spatialScalesSanityCheck(params_TYP * params, FS_TYP * FS)
{
	MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // =================
	if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << endl << "* * * * * * * * * * * * CHECKING VALIDITY OF HYBRID MODEL FOR THE SIMULATED PLASMA * * * * * * * * * * * * * * * * * *" << endl;
        cout << "Electron skin depth to grid size ratio: " << scientific << FS->electronSkinDepth/params->mesh.DX << fixed << endl;
        cout << "* * * * * * * * * * * * VALIDITY OF HYBRID MODEL FOR THE SIMULATED PLASMA CHECKED  * * * * * * * * * * * * * * * * * *" << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Check that DX is larger than the electron skin depth, otherwise, abort simulation:
        // ==================================================================================
	if (params->mesh.DX <= FS->electronSkinDepth)
    {
        cout << "ERROR: Grid size violates assumptions of hybrid model for the plasma -- lenght scales smaller than the electron skind depth can not be resolved." << endl;
        cout << "ABORTING SIMULATION..." << endl;

        MPI_Abort(MPI_COMM_WORLD,-103);
	}
}

void units_TYP::defineTimeStep(params_TYP * params, vector<ionSpecies_TYP> * IONS)
{
	/*
    CFL condition for Whistler waves has not been included as we are currently
    not solving the magnetic field via Ampere's law.
    */

    // Print to terminal:
    // ==================
	if (params->mpi.IS_PARTICLES_ROOT)
    {
        cout << endl << "* * * * * * * * * * * * COMPUTING SIMULATION TIME STEP * * * * * * * * * * * * * * * * * *" << endl;
    }

    // Define variables:
    // =================
	double ionsMaxVel(0.0);	// Maximum speed of simulated ions
	double DT(0.0); 	// Time step defined by user
	double DT_CFL_I(0.0);	// Minimum time step defined by CFL condition for ions
	bool CFL_I(false);

	// Time step given by user:
    // ========================
	DT = params->DTc*params->ionGyroPeriod;

	// CFL condition for ions:
    // =======================
	if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
        {
            for (int ss=0; ss<params->numberOfParticleSpecies; ss++)
            {
                    vec V = sqrt( pow(IONS->at(ss).V_p.col(0), 2.0) + pow(IONS->at(ss).V_p.col(1), 2.0) );

                    ionsMaxVel = (ionsMaxVel < V.max()) ? V.max() : ionsMaxVel;
            }

            // Minimum time step required by CFL condition for ions:
            DT_CFL_I = params->mesh.DX/ionsMaxVel;

            // We gather DT_CFL_I from all MPI processes in particles communicator:
            double * DT_CFL_I_MPI;

            DT_CFL_I_MPI = (double*)malloc( params->mpi.MPIS_PARTICLES*sizeof(double) );

            MPI_Allgather(&DT_CFL_I, 1, MPI_DOUBLE, DT_CFL_I_MPI, 1, MPI_DOUBLE, params->mpi.COMM);

            // Sort the DT_CFL_I_MPI values:
            // =============================
            if (params->mpi.IS_PARTICLES_ROOT)
            {
                // Find maximum value of DT_CFL_I accross all MPIs:
                for (int ii=0; ii<params->mpi.MPIS_PARTICLES; ii++)
                {
                    DT_CFL_I = (DT_CFL_I > *(DT_CFL_I_MPI + ii)) ? *(DT_CFL_I_MPI + ii) : DT_CFL_I;
                }

                // Check against user defined time step:
                if (DT > DT_CFL_I)
                {
                        DT = DT_CFL_I;
                        CFL_I = true;
                }

				//Assign final DT for the simulation
				//==================================
				params->DT = DT;

                params->timeIterations = (int)ceil( params->simulationTime*params->ionGyroPeriod/params->DT );

                params->outputCadenceIterations = (int)ceil( params->outputCadence*params->ionGyroPeriod/params->DT );
            }

            // Deallocate memory:
            free(DT_CFL_I_MPI);
	}

	// Broadcast correct time step, time iterations in simulation, and cadence for generating outputs:
    // ===============================================================================================
	MPI_Bcast(&params->DT, 1, MPI_DOUBLE, params->mpi.PARTICLES_ROOT_WORLD_RANK, MPI_COMM_WORLD);

	MPI_Bcast(&params->timeIterations, 1, MPI_INT, params->mpi.PARTICLES_ROOT_WORLD_RANK, MPI_COMM_WORLD);

	MPI_Bcast(&params->outputCadenceIterations, 1, MPI_INT, params->mpi.PARTICLES_ROOT_WORLD_RANK, MPI_COMM_WORLD);

	// Print to terminal a brief summary:
    // ==================================
	if (params->mpi.IS_PARTICLES_ROOT)
        {
            if (CFL_I)
            {
                cout << "+ Simulation time step defined by CFL condition for IONS" << endl;
            }
            else
            {
                cout << "+ Simulation time step defined by USER: " <<  scientific << DT << fixed  << endl;
            }

            cout << "+ Time step defined by CFL condition for ions: " << scientific << DT_CFL_I << fixed << endl;
            cout << "+ Time step used in simulation: " << scientific << params->DT << fixed << endl;
            cout << "+ Time steps in simulation: " << params->timeIterations << endl;
            cout << "+ Simulation time: " << scientific << params->DT*params->timeIterations << fixed << " s" << endl;
            cout << "+ Simulation time: " << params->DT*params->timeIterations/params->ionGyroPeriod << " gyroperiods" << endl;
            cout << "+ Cadence for saving outputs: " << params->outputCadenceIterations << endl;
            cout << "+ Number of outputs: " << floor(params->timeIterations/params->outputCadenceIterations) + 1 << endl;
            cout << "+ Cadence for checking stability: " << params->rateOfChecking << endl;
            cout << "* * * * * * * * * * * * * * * TIME STEP COMPUTED * * * * * * * * * * * * * * * * * * * * *" << endl;
	}
}

void units_TYP::normalizeVariables(params_TYP * params, vector<ionSpecies_TYP> * IONS, fields_TYP * fields, const CS_TYP * CS)
{
	// Normalizing physical constants:
    // =========================================================================
	F_E_DS /= CS->charge; 					// Dimensionless electron charge
	F_ME_DS /= CS->mass; 					// Dimensionless electron charge
	F_MU_DS /= CS->vacuumPermeability; 		// Dimensionless vacuum permeability
	F_EPSILON_DS /= CS->vacuumPermittivity;	// Dimensionless vacuum permittivity
	F_C_DS /= CS->velocity; 				// Dimensionless speed of light

	// Normalizing "params":
    // =========================================================================
	// time increment:
	// ---------------
	params->DT /= CS->time;

    /*
    We need to understand how ne is normalized and how A_0 needs
    to be normalized too. Ideally we would like ne to be in m^-3
    and A_0 in m^2

    Consider that if we normalize A_0, we will need to add its contribution once
    Save the output data, in other words, multiply length*area
    */

	// Characteristic values:
	// ----------------------
	params->CV.ne /= CS->density;
	params->CV.Te /= CS->temperature;
	params->CV.B /= CS->bField;
	params->CV.Tpar /= CS->temperature;
	params->CV.Tper /= CS->temperature;

	// Fluid initial conditions:
	// -------------------------
	params->f_IC.ne /= CS->density;
	params->f_IC.Te /= CS->temperature;

	// Electromagnetic fields initial conditions:
	// -----------------------------------------
	params->em_IC.BX 	 	 /= CS->bField;
	params->em_IC.BY     	 /= CS->bField;
	params->em_IC.BZ         /= CS->bField;
	params->em_IC.Bx_profile /= CS->bField;

	params->em_IC.EX 	 	 /= CS->eField;
	params->em_IC.EY     	 /= CS->eField;
	params->em_IC.EZ         /= CS->eField;
	params->em_IC.Ex_profile /= CS->eField;

	// Geometry:
	// ---------
	params->geometry.r1 /= CS->length;
	params->geometry.r2 /= CS->length;
    // params->geometry.A_0 /= pow(CS->length,2);

	// Fundamental scales:
	// -------------------
	params->ionLarmorRadius /= CS->length;
	params->ionSkinDepth /= CS->length;
	params->ionGyroPeriod /= CS->time;

	// Mesh quantities:
	// ----------------
	params->mesh.nodesX = params->mesh.nodesX/CS->length;
	params->mesh.nodesY = params->mesh.nodesY/CS->length;
	params->mesh.nodesZ = params->mesh.nodesZ/CS->length;
	params->mesh.DX /= CS->length;
	params->mesh.DY /= CS->length;
	params->mesh.DZ /= CS->length;
	params->mesh.LX /= CS->length;
	params->mesh.LY /= CS->length;
	params->mesh.LZ /= CS->length;

	// Normalizing IONS:
    // =========================================================================
	for(int ii=0;ii<IONS->size();ii++)
    {
		// Particle boundary conditions:
		// ----------------------------
		IONS->at(ii).p_BC.mean_x  /= CS->length;
		IONS->at(ii).p_BC.sigma_x /= CS->length;
		IONS->at(ii).p_BC.T /= CS->temperature;
		IONS->at(ii).p_BC.E /= CS->temperature;
        IONS->at(ii).p_BC.G *= CS->time;

		// Particle parameters:
		// -------------------
		IONS->at(ii).Q /= CS->charge;
		IONS->at(ii).M /= CS->mass;

		// Particle initial conditions:
		// ----------------------------
		IONS->at(ii).p_IC.Tpar /= CS->temperature;
		IONS->at(ii).p_IC.Tper /= CS->temperature;
		IONS->at(ii).p_IC.Tper_profile /= CS->temperature;
		IONS->at(ii).p_IC.Tpar_profile /= CS->temperature;
		IONS->at(ii).p_IC.densityFraction_profile /= CS->density;

        // Scales:
		IONS->at(ii).LarmorRadius /= CS->length;
		IONS->at(ii).VTpar /= CS->velocity;
		IONS->at(ii).VTper /= CS->velocity;
		IONS->at(ii).Wc *= CS->time;
		IONS->at(ii).Wp *= CS->time;

		if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    	{
                    IONS->at(ii).X_p = IONS->at(ii).X_p/CS->length;
                    IONS->at(ii).V_p = IONS->at(ii).V_p/CS->velocity;
		}
	}

    // Normalizing "fields":
    // =========================================================================
	fields->EX_m   /=  CS->eField;
	fields->BX_m   /=  CS->bField;
    fields->dBX_m  /= (CS->bField/pow(CS->length,1));
    fields->ddBX_m /= (CS->bField/pow(CS->length,2));
}
