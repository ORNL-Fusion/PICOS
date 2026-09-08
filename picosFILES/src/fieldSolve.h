#ifndef H_FIELDS_SOLVER
#define H_FIELDS_SOLVER

#include <iostream>
#include <cmath>
#include <vector>

#define ARMA_ALLOW_FAKE_GCC
#include "armadillo"
#include "types.h"

#include "mpi_main.h"

using namespace std;
using namespace arma;

class fields_solver_TYP
{

  	int NX_S; // Number of grid cells per subdomain, including 2 ghost cells.
  	int NX_T; // Number of grid cells in entire simulation domain, including 2 ghost cells.
  	int NX_R; // Number of grid cells in entire simulation domain, not including ghost cells.

  	// Electron density:
  	arma::vec ne;    // Current time tt
  	arma::vec ne_;   // tt - 1
  	arma::vec ne__;  // tt - 2
  	arma::vec ne___; // tt - 3

    // Electron temperature:
    arma::vec Te;

    // Electron pressure:
    arma::vec Pe;

    // Electron pressure gradient:
    arma::vec dPe;

  	// Electron density gradient:
  	//arma::vec dne;

    // Electric field:
    arma::vec EX_m;

    // Electrostatic Poisson fields:
    arma::vec Phi_m;
    arma::vec chargeDensity;

    // Reformulated Poisson fields:
    arma::vec ionDensity;
    arma::vec electronDensity;
    arma::vec stressDifference;
    arma::vec divStressDifference;
    arma::vec plasmaFrequencySquared;

  	// Grid cell increment
  	double dx;

  	// MPI functions:
  	void MPI_Allgathervec(const params_TYP * params, arma::vec * field);

  	void MPI_SendVec(const params_TYP * params, arma::vec * v);

    // Ghost cells:
    void fillGhosts(arma::vec * C);

    void fillPeriodicGhosts(arma::vec * C);

    void fill4Ghosts(arma::vec * v);

    // Smoothing:
    void smooth(arma::vec * v, double as);

    void smoothPeriodic(arma::vec * v, double as);

    void advanceEfieldOhmLaw(const params_TYP * params, fields_TYP * fields, CS_TYP * CS, vector<ionSpecies_TYP> * IONS, electrons_TYP * electrons);

    void advanceEfieldPoisson(const params_TYP * params, fields_TYP * fields, CS_TYP * CS, vector<ionSpecies_TYP> * IONS);

    void advanceEfieldReformulatedPoisson(const params_TYP * params, fields_TYP * fields, CS_TYP * CS, vector<ionSpecies_TYP> * IONS);

    void solveDirichletPoisson(const params_TYP * params, const arma::vec * rho, arma::vec * phi) const;

    void solvePeriodicPoisson(const params_TYP * params, const arma::vec * rho, arma::vec * phi) const;

    double poissonBoundaryPotential(const params_TYP * params, bool rightBoundary) const;

	public:

  	fields_solver_TYP(){};

  	fields_solver_TYP(const params_TYP * params, CS_TYP * CS);

  	//void advanceBField(const params_TYP * params, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

  	void advanceEfield(const params_TYP * params, fields_TYP * fields, CS_TYP * CS, vector<ionSpecies_TYP> * IONS, electrons_TYP * electrons);
};

#endif
