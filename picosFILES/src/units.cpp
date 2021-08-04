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
