#ifndef H_OUTPUTHDF5
#define H_OUTPUTHDF5

#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif

#include <string>
#include <cmath>

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>
#include "types.h"

#include "H5Cpp.h"

using namespace H5;
using namespace std;
using namespace arma;

class HDF_TYP
{

	#ifdef HDF5_DOUBLE
		#define HDF_TYPE PredType::NATIVE_DOUBLE
		#define CPP_TYPE double
	#elif defined HDF5_FLOAT
		#define HDF_TYPE PredType::NATIVE_FLOAT
		#define CPP_TYPE float
	#endif

	/*
	void MPI_Allgathervec(const params_TYP * params, arma::vec * field);

	void MPI_Allgathermat(const params_TYP * params, arma::mat * field);
	*/

	void saveToHDF5(H5File * file, string name, int * value);

	void saveToHDF5(H5File * file, string name, CPP_TYPE * value);

	void saveToHDF5(Group * group, string name, int * value);

	void saveToHDF5(Group * group, string name, CPP_TYPE * value);

	void saveToHDF5(H5File * file, string name, std::vector<int> * values);

	void saveToHDF5(H5File * file, string name, std::vector<CPP_TYPE> * values);

	void saveToHDF5(H5File * file, string name, arma::ivec * values);

	void saveToHDF5(Group * group, string name, arma::ivec * values);

	void saveToHDF5(H5File * file, string name, arma::vec * values);

	void saveToHDF5(Group * group, string name, arma::vec * values);

	void saveToHDF5(Group * group, string name, arma::fvec * values);

	void saveToHDF5(H5File * file, string name, arma::imat * values);

	void saveToHDF5(Group * group, string name, arma::imat * values);

	void saveToHDF5(Group * group, string name, arma::mat * values);

	void saveToHDF5(Group * group, string name, arma::fmat * values);


	void saveIonsVariables(const params_TYP * params, const vector<ionSpecies_TYP> * IONS, const CS_TYP * CS, const Group * group_iteration);

	void saveFieldsVariables(const params_TYP * params, fields_TYP * fields, const CS_TYP * CS, const Group * group_iteration);

public:

	HDF_TYP(params_TYP * params, FS_TYP * FS, vector<ionSpecies_TYP> * IONS);

	void saveOutputs(const params_TYP * params, const vector<ionSpecies_TYP> * IONS, fields_TYP * fields, const CS_TYP * CS, const int it, double totalTime);
};

#endif
