#ifndef H_INITIALIZE
#define H_INITIALIZE

// Intrinsic header files:
// =======================
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <typeinfo>

// Armadillo header:
// =================
#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

// User-defined headers:
// =====================
#include "types.h"
//#include "quietStart.h"
//#include "PIC.h"
#include "mpi_main.h"

using namespace std;
using namespace arma;

class init_TYP
{
	vector<string> split(const string& str, const string& delim);

	map<string, string> readTextFile(string *  inputFile);

	void allocateParticleDefinedIonArrays(const params_TYP * params, ionSpecies_TYP * IONS);

	void allocateMeshDefinedIonArrays(const params_TYP * params, ionSpecies_TYP * IONS);

public:

	init_TYP(params_TYP * params, int argc, char* argv[]);

	void calculateMeshParams(params_TYP * params);

	void readInputFile(params_TYP * params);

  	void readIonPropertiesFile(params_TYP * params, vector<ionSpecies_TYP> * IONS);

	void calculateDerivedQuantities(params_TYP * params, vector<ionSpecies_TYP> * IONS);

  	void readInitialConditionProfiles(params_TYP * params, vector<ionSpecies_TYP> * IONS);

	void setupIonsInitialCondition(const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

	void initializeFields(params_TYP * params, fields_TYP * fields);

	void allocateMemoryIons(params_TYP * params, vector<ionSpecies_TYP> * IONS);

};

#endif
