#ifndef H_RFOPERATOR
#define H_RFOPERATOR

#include <iostream>
#include <cmath>
#include <vector>

#include "armadillo"
#include "types.h"
#include "mpi_main.h"
#include "omp.h"

using namespace std;
using namespace arma;

class RF_Operator_TYP
{
    void calculateResNum_AllSpecies( params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

    void calculateResNum(int ii,  params_TYP * params, CS_TYP * CS, fields_TYP * fields, ionSpecies_TYP * IONS);

    void checkResNumAndFlag_AllSpecies( params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

    void calculateRfTerms_AllSpecies( params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

    void calculateRfTerms(int ii,  params_TYP * params, CS_TYP * CS, fields_TYP * fields, ionSpecies_TYP * IONS);

    void calculatePowerPerUnitErf_AllSpecies( params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

    void calculateAbsorbedPower_AllSpecies(params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

    void ApplyRfOperator_AllSpecies(params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

    void MPI_AllreduceDouble( params_TYP * params, double * v);

  public:

    RF_Operator_TYP( params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

    void ApplyRfHeating_AllSpecies( params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);
};

#endif
