#ifndef H_PIC
#define H_PIC

// System libraries:
#include <iostream>
#include <cmath>
#include <vector>

/*
#ifdef __linux__
#include <chrono>//C++11 standard
#include <random>//C++11 standard

#elif __APPLE__
//define something else for apple compiler
#endif
*/

// User-defined libraries:
#include "armadillo"
#include "types.h"

// Parallelization libraries:
#include <omp.h>
#include "mpi_main.h"

/*
//#include <cfenv>
//#pragma STDC FENV_ACCESS ON
*/

using namespace std;
using namespace arma;

class PIC_TYP
{
  /*
  int NX;
	int NY;
	int NZ;
	int ix;
	int iy;
	int iz;
	arma::vec x;
	uvec logic;
  */

protected:

  // MPI methods:
	//void MPI_AllreduceVec(const params_TYP * params, arma::vec * v);

	void MPI_SendVec(const params_TYP * params, arma::vec * v);

	void MPI_ReduceVec(const params_TYP * params, arma::vec * v);

	void MPI_Allgathervec(const params_TYP * params, arma::vec * field);

	void MPI_Recvvec(const params_TYP * params, arma::vec * field);

  // Ghost contributions:
	void include4GhostsContributions(arma::vec * v);

  void fillGhosts(arma::vec * C);

	void fill4Ghosts(arma::vec * v);

  // Smoothing:
	void smooth(arma::vec * v, double as);

  void interpolateScalarField(const params_TYP * params, ionSpecies_TYP * IONS, arma::vec * F_m, arma::vec * F_p);

  void interpEM(const params_TYP * params, CS_TYP * CS, const ionSpecies_TYP * IONS, const fields_TYP * fields, arma::rowvec * ZN, arma::rowvec * EM);

  void calculateF(const params_TYP * params, CS_TYP * CS, const ionSpecies_TYP * IONS, arma::rowvec * ZN, arma::rowvec * EM, arma::rowvec * F);

// void interpolateVectorField(const params_TYP * params, const ionSpecies_TYP * IONS, vfield_vec * field, arma::mat * F);

	// void interpolateElectromagneticFields(const params_TYP * params, const ionSpecies_TYP * IONS, fields_TYP * fields, arma::mat * E, arma::mat * B);

	//void interpolateVectorField(const params_TYP * params, const twoDimensional::ionSpecies * IONS, vfield_mat * field, arma::mat * F);

	//void interpolateElectromagneticFields(const params_TYP * params, const twoDimensional::ionSpecies * IONS, twoDimensional::fields * EB, arma::mat * E, arma::mat * B);

	void eim(const params_TYP * params, fields_TYP * fields, ionSpecies_TYP * IONS);

	void calculateIonMoments(const params_TYP * params, fields_TYP * fields, ionSpecies_TYP * IONS);

	void calculateDerivedIonMoments(const params_TYP * params, ionSpecies_TYP * IONS);

  public:

	PIC_TYP();

	void assignCell(const params_TYP * params, fields_TYP * fields, ionSpecies_TYP * IONS);

  void advanceParticles(const params_TYP * params, CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

  void extrapolateIonsMoments(const params_TYP * params, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

  /*
	void advanceIonsVelocity(const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS, const double DT);

	void advanceIonsPosition(const params_TYP * params,  fields_TYP * fields, vector<ionSpecies_TYP> * IONS, const double DT);
  */


};

#endif
