#ifndef H_UNITS
#define H_UNITS

#include <iostream>
#include <vector>
#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>
#include <cmath>

#include "types.h"
#include "mpi_main.h"

using namespace std;
using namespace arma;

class units_TYP
{

	void broadcastCharacteristicScales(params_TYP * params, CS_TYP * CS);

	void broadcastFundamentalScales(params_TYP * params, FS_TYP * FS);

public:

	units_TYP(){};

	void spatialScalesSanityCheck(params_TYP * params, FS_TYP * FS);

	void defineTimeStep(params_TYP * params, vector<ionSpecies_TYP> * IONS);

	void calculateFundamentalScales(params_TYP * params, vector<ionSpecies_TYP> * IONS, FS_TYP * FS);

	void defineCharacteristicScales(params_TYP * params, vector<ionSpecies_TYP> * IONS, CS_TYP * CS);

	void normalizeVariables(params_TYP * params, vector<ionSpecies_TYP> * IONS, fields_TYP * fields, const CS_TYP * CS);

	void defineCharacteristicScalesAndBcast(params_TYP * params, vector<ionSpecies_TYP> * IONS, CS_TYP * CS);

	void calculateFundamentalScalesAndBcast(params_TYP * params, vector<ionSpecies_TYP> * IONS, FS_TYP * FS);

};

#endif
