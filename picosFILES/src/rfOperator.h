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
    void calculateResNum_AllSpecies(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

    void calculateResNum(int ii, const params_TYP * params, CS_TYP * CS, fields_TYP * fields, ionSpecies_TYP * IONS);

    void checkResNumAndFlag_AllSpecies(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

    void calculateRfTerms_AllSpecies(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

    void calculateRfTerms(int ii, const params_TYP * params, CS_TYP * CS, fields_TYP * fields, ionSpecies_TYP * IONS);

  public:

    RF_Operator_TYP(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

    void ApplyRfHeating(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);
};

#endif
