#ifndef H_PARTICLEBOUNDARYCONDITIONS
#define H_PARTICLEBOUNDARYCONDITIONS

#include <iostream>
#include <cmath>
#include <vector>

#define ARMA_ALLOW_FAKE_GCC
#include "armadillo"
#include "types.h"
#include "mpi_main.h"
#include "omp.h"

using namespace std;
using namespace arma;

class particleBC_TYP
{
private:

    //void particleReinjection(int ii, const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, ionSpecies_TYP * IONS);

    void particleReinjection(int ii, const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, ionSpecies_TYP * IONS);

    void MPI_AllreduceDouble(const params_TYP * params, double * v);

    template <typename vec_TYP> void MPI_OMP_AllreduceVec(const params_TYP * params, vec_TYP * V, double * S);

    void checkBoundaryAndFlag(const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

    void calculateParticleWeight(const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

public:
    particleBC_TYP();

    void applyParticleReinjection(const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);
};

#endif
