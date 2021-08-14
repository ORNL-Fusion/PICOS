#include "outputHDF5.h"

void HDF_TYP::MPI_Allgathervec(const params_TYP * params, arma::vec * field)
{
	unsigned int iIndex = params->mpi.iIndex;
	unsigned int fIndex = params->mpi.fIndex;

	arma::vec recvbuf(params->mesh.NX_IN_SIM);
	arma::vec sendbuf(params->mesh.NX_PER_MPI);

	//Allgather for x-component
	sendbuf = field->subvec(iIndex, fIndex);
	MPI_Allgather(sendbuf.memptr(), params->mesh.NX_PER_MPI, MPI_DOUBLE, recvbuf.memptr(), params->mesh.NX_PER_MPI, MPI_DOUBLE, params->mpi.MPI_TOPO);
	field->subvec(1, params->mesh.NX_IN_SIM) = recvbuf;
}

void HDF_TYP::MPI_Allgathermat(const params_TYP * params, arma::mat * field)
{
	unsigned int irow = params->mpi.irow;
	unsigned int frow = params->mpi.frow;

	unsigned int icol = params->mpi.icol;
	unsigned int fcol = params->mpi.fcol;

	arma::vec recvbuf = zeros(params->mesh.NX_IN_SIM*params->mesh.NY_IN_SIM);
	arma::vec sendbuf = zeros(params->mesh.NX_PER_MPI*params->mesh.NY_PER_MPI);

	//Allgather for x-component
	sendbuf = vectorise(field->submat(irow,icol,frow,fcol));
	MPI_Allgather(sendbuf.memptr(), params->mesh.NUM_CELLS_PER_MPI, MPI_DOUBLE, recvbuf.memptr(), params->mesh.NUM_CELLS_PER_MPI, MPI_DOUBLE, params->mpi.MPI_TOPO);

	for (int mpis=0; mpis<params->mpi.NUMBER_MPI_DOMAINS; mpis++)
  {
		unsigned int ie = params->mesh.NX_PER_MPI*params->mesh.NY_PER_MPI*mpis;
		unsigned int fe = params->mesh.NX_PER_MPI*params->mesh.NY_PER_MPI*(mpis+1) - 1;

		unsigned int ir = *(params->mpi.MPI_CART_COORDS.at(mpis))*params->mesh.NX_PER_MPI + 1;
		unsigned int fr = ( *(params->mpi.MPI_CART_COORDS.at(mpis)) + 1)*params->mesh.NX_PER_MPI;

		unsigned int ic = *(params->mpi.MPI_CART_COORDS.at(mpis)+1)*params->mesh.NY_PER_MPI + 1;
		unsigned int fc = ( *(params->mpi.MPI_CART_COORDS.at(mpis)+1) + 1)*params->mesh.NY_PER_MPI;

		field->submat(ir,ic,fr,fc) = reshape(recvbuf.subvec(ie,fe), params->mesh.NX_PER_MPI, params->mesh.NY_PER_MPI);
	}
}


// Function to save a single integer value
void HDF_TYP::saveToHDF5(H5File * file, string name, int * value)
{
	H5std_string nameSpace( name );
	int data[1] = {*value};
	hsize_t dims[1] = {1};
	DataSpace * dataspace = new DataSpace(1, dims);
	DataSet * dataset = new DataSet(file->createDataSet( nameSpace, PredType::NATIVE_INT, *dataspace ));

	dataset->write( data, PredType::NATIVE_INT);

	delete dataspace;
	delete dataset;
}


// Function to save a single CPP_TYPE value
void HDF_TYP::saveToHDF5(H5File * file, string name, CPP_TYPE * value)
{
	H5std_string nameSpace( name );
	CPP_TYPE data[1] = {*value};
	hsize_t dims[1] = {1};
	DataSpace * dataspace = new DataSpace(1, dims);
	DataSet * dataset = new DataSet(file->createDataSet( nameSpace, HDF_TYPE, *dataspace ));

	dataset->write( data, HDF_TYPE);

	delete dataspace;
	delete dataset;
}


// Function to save a single integer value (to a HDF5 group)
void HDF_TYP::saveToHDF5(Group * group, string name, int * value)
{
	H5std_string nameSpace( name );
	int data[1] = {*value};
	hsize_t dims[1] = {1};
	DataSpace * dataspace = new DataSpace(1, dims);
	DataSet * dataset = new DataSet(group->createDataSet( nameSpace, PredType::NATIVE_INT, *dataspace ));

	dataset->write( data, PredType::NATIVE_INT);

	delete dataspace;
	delete dataset;
}


// Function to save a single CPP_TYPE value (to a HDF5 group)
void HDF_TYP::saveToHDF5(Group * group, string name, CPP_TYPE * value)
{
	H5std_string nameSpace( name );
	CPP_TYPE data[1] = {*value};
	hsize_t dims[1] = {1};
	DataSpace * dataspace = new DataSpace(1, dims);
	DataSet * dataset = new DataSet(group->createDataSet( nameSpace, HDF_TYPE, *dataspace ));

	dataset->write( data, HDF_TYPE);

	delete dataspace;
	delete dataset;
}


// Function to save a vector of int values
void HDF_TYP::saveToHDF5(H5File * file, string name, std::vector<int> * values)
{
	H5std_string nameSpace( name );
	unsigned long long int size = (unsigned long long int)values->size();

	int * data;
   	data = new int[size];
        std::copy(values->begin(), values->end(), data);

	hsize_t dims[1] = {size};
	DataSpace * dataspace = new DataSpace(1, dims);
	DataSet * dataset = new DataSet(file->createDataSet( nameSpace, PredType::NATIVE_INT, *dataspace ));

	dataset->write( data, PredType::NATIVE_INT);

	delete dataspace;
	delete dataset;
}


// Function to save a vector of CPP_TYPE values
void HDF_TYP::saveToHDF5(H5File * file, string name, std::vector<CPP_TYPE> * values)
{
	H5std_string nameSpace( name );
	unsigned long long int size = (unsigned long long int)values->size();

	CPP_TYPE * data;
   	data = new CPP_TYPE[size];
    std::copy(values->begin(), values->end(), data);

	hsize_t dims[1] = {size};
	DataSpace * dataspace = new DataSpace(1, dims);
	DataSet * dataset = new DataSet(file->createDataSet( nameSpace, HDF_TYPE, *dataspace ));

	dataset->write( data, HDF_TYPE);

	delete dataspace;
	delete dataset;
}


// Function to save an Armadillo ivec vector (to a HDF5 file)
void HDF_TYP::saveToHDF5(H5File * file, string name, arma::ivec * values)
{
	H5std_string nameSpace( name );

	hsize_t dims[1] = {(hsize_t)values->n_elem};
	DataSpace * dataspace = new DataSpace(1, dims);
	DataSet * dataset = new DataSet(file->createDataSet( nameSpace, PredType::NATIVE_INT, *dataspace ));

	dataset->write( values->memptr(), PredType::NATIVE_INT);

	delete dataspace;
	delete dataset;
}


// Function to save an Armadillo ivec vector (to a HDF5 group)
void HDF_TYP::saveToHDF5(Group * group, string name, arma::ivec * values)
{
	H5std_string nameSpace( name );

	hsize_t dims[1] = {(hsize_t)values->n_elem};
	DataSpace * dataspace = new DataSpace(1, dims);
	DataSet * dataset = new DataSet(group->createDataSet( nameSpace, PredType::NATIVE_INT, *dataspace ));

	dataset->write( values->memptr(), PredType::NATIVE_INT);

	delete dataspace;
	delete dataset;
}


// Function to save an Armadillo vec vector
void HDF_TYP::saveToHDF5(H5File * file, string name, arma::vec * values)
{
	H5std_string nameSpace( name );

	hsize_t dims[1] = {(hsize_t)values->n_elem};
	DataSpace * dataspace = new DataSpace(1, dims);
	DataSet * dataset = new DataSet(file->createDataSet( nameSpace, HDF_TYPE, *dataspace ));

	dataset->write( values->memptr(), HDF_TYPE);

	delete dataspace;
	delete dataset;
}


// Function to save an Armadillo vec vector (to a HDF5 group)
void HDF_TYP::saveToHDF5(Group * group, string name, arma::vec * values)
{
	H5std_string nameSpace( name );

	hsize_t dims[1] = {(hsize_t)values->n_elem};
	DataSpace * dataspace = new DataSpace(1, dims);
	DataSet * dataset = new DataSet(group->createDataSet( nameSpace, HDF_TYPE, *dataspace ));

	dataset->write( values->memptr(), HDF_TYPE);

	delete dataspace;
	delete dataset;
}


// Function to save an Armadillo fvec vector (to a HDF5 group)
void HDF_TYP::saveToHDF5(Group * group, string name, arma::fvec * values)
{
	H5std_string nameSpace( name );

	hsize_t dims[1] = {(hsize_t)values->n_elem};
	DataSpace * dataspace = new DataSpace(1, dims);
	DataSet * dataset = new DataSet(group->createDataSet( nameSpace, HDF_TYPE, *dataspace ));

	dataset->write( values->memptr(), HDF_TYPE);

	delete dataspace;
	delete dataset;
}


// Function to save an Armadillo imat vector (to a HDF5 file)
void HDF_TYP::saveToHDF5(H5File * file, string name, arma::imat * values)
{
	H5std_string nameSpace( name );

	hsize_t dims[2] = {(hsize_t)values->n_cols, (hsize_t)values->n_rows};
	DataSpace * dataspace = new DataSpace(2, dims);
	DataSet * dataset = new DataSet(file->createDataSet( nameSpace, PredType::NATIVE_INT, *dataspace ));

	dataset->write( values->memptr(), PredType::NATIVE_INT);

	delete dataspace;
	delete dataset;
}


// Function to save an Armadillo imat vector (to a HDF5 group)
void HDF_TYP::saveToHDF5(Group * group, string name, arma::imat * values)
{
	H5std_string nameSpace( name );

	hsize_t dims[2] = {(hsize_t)values->n_cols, (hsize_t)values->n_rows};
	DataSpace * dataspace = new DataSpace(2, dims);
	DataSet * dataset = new DataSet(group->createDataSet( nameSpace, PredType::NATIVE_INT, *dataspace ));

	dataset->write( values->memptr(), PredType::NATIVE_INT);

	delete dataspace;
	delete dataset;
}


// Function to save an Armadillo vec vector (to a HDF5 group)
void HDF_TYP::saveToHDF5(Group * group, string name, arma::mat * values)
{
	H5std_string nameSpace( name );

	hsize_t dims[2] = {(hsize_t)values->n_cols, (hsize_t)values->n_rows};
	DataSpace * dataspace = new DataSpace(2, dims);
	DataSet * dataset = new DataSet(group->createDataSet( nameSpace, HDF_TYPE, *dataspace ));

	dataset->write( values->memptr(), HDF_TYPE);

	delete dataspace;
	delete dataset;
}


// Function to save an Armadillo vec vector (to a HDF5 group)
void HDF_TYP::saveToHDF5(Group * group, string name, arma::fmat * values)
{
	H5std_string nameSpace( name );

	hsize_t dims[2] = {(hsize_t)values->n_cols, (hsize_t)values->n_rows};
	DataSpace * dataspace = new DataSpace(2, dims);
	DataSet * dataset = new DataSet(group->createDataSet( nameSpace, HDF_TYPE, *dataspace ));

	dataset->write( values->memptr(), HDF_TYPE);

	delete dataspace;
	delete dataset;
}


// Constructor:
// =============================================================================
HDF_TYP::HDF_TYP(params_TYP * params, FS_TYP * FS, vector<ionSpecies_TYP> * IONS)
{
    if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
      try
      {
        string name;
        string path;

        int int_value;
        CPP_TYPE cpp_type_value;
        std::vector<CPP_TYPE> vector_values;

        arma::vec vec_values;
        arma::fvec fvec_values;

        path = params->PATH + "/HDF5/";
        name = path + "main.h5";
        const H5std_string	FILE_NAME( name );
        name.clear();

        Exception::dontPrint();

        // Create a new file using the default property lists.
        H5File * outputFile = new H5File( FILE_NAME, H5F_ACC_TRUNC );

        /*
        name = "dimensionality";
        saveToHDF5(outputFile, name, &params->dimensionality);
        name.clear();
        */

        name = "smoothingParameter";
        cpp_type_value = params->smoothingParameter;
        saveToHDF5(outputFile, name, &cpp_type_value);
        name.clear();

        name = "numberOfRKIterations";
        saveToHDF5(outputFile, name, &params->numberOfRKIterations);
        name.clear();

        name = "filtersPerIterationFields";
        saveToHDF5(outputFile, name, &params->filtersPerIterationFields);
        name.clear();

        name = "numOfMPIs";
        saveToHDF5(outputFile, name, &params->mpi.NUMBER_MPI_DOMAINS);
        name.clear();

        name = "numMPIsParticles";
        saveToHDF5(outputFile, name, &params->mpi.MPIS_PARTICLES);
        name.clear();


        name = "numMPIsFields";
        saveToHDF5(outputFile, name, &params->mpi.MPIS_FIELDS);
        name.clear();

        // Fundamental scales group
        Group * group_scales = new Group( outputFile->createGroup( "/scales" ) );

        name = "electronSkinDepth";
        cpp_type_value = FS->electronSkinDepth;
        saveToHDF5(group_scales, name, &cpp_type_value);
        name.clear();

        name = "electronGyroPeriod";
        cpp_type_value = FS->electronGyroPeriod;
        saveToHDF5(group_scales, name, &cpp_type_value);
        name.clear();

        name = "electronGyroRadius";
        cpp_type_value = FS->electronGyroRadius;
        saveToHDF5(group_scales, name, &cpp_type_value);
        name.clear();

        name = "ionGyroRadius";
        #ifdef HDF5_DOUBLE
        vec_values = zeros(params->numberOfParticleSpecies);
        for (int ss=0; ss<params->numberOfParticleSpecies; ss++)
                vec_values(ss) = FS->ionGyroRadius[ss];
        saveToHDF5(group_scales, name, &vec_values);
        #elif defined HDF5_FLOAT
        fvec_values = zeros<fvec>(params->numberOfParticleSpecies);
        for (int ss=0; ss<params->numberOfParticleSpecies; ss++)
                fvec_values(ss) = (float)FS->ionGyroRadius[ss];
        saveToHDF5(group_scales, name, &fvec_values);
        #endif
        name.clear();

        name = "ionGyroPeriod";
        #ifdef HDF5_DOUBLE
        vec_values = zeros(params->numberOfParticleSpecies);
        for (int ss=0; ss<params->numberOfParticleSpecies; ss++)
                vec_values(ss) = FS->ionGyroPeriod[ss];
        saveToHDF5(group_scales, name, &vec_values);
        #elif defined HDF5_FLOAT
        fvec_values = zeros<fvec>(params->numberOfParticleSpecies);
        for (int ss=0; ss<params->numberOfParticleSpecies; ss++)
                fvec_values(ss) = (float)FS->ionGyroPeriod[ss];
        saveToHDF5(group_scales, name, &fvec_values);
        #endif
        name.clear();

        name = "ionSkinDepth";
        #ifdef HDF5_DOUBLE
        vec_values = zeros(params->numberOfParticleSpecies);
        for (int ss=0; ss<params->numberOfParticleSpecies; ss++)
                vec_values(ss) = FS->ionSkinDepth[ss];
        saveToHDF5(group_scales, name, &vec_values);
        #elif defined HDF5_FLOAT
        fvec_values = zeros<fvec>(params->numberOfParticleSpecies);
        for (int ss=0; ss<params->numberOfParticleSpecies; ss++)
                fvec_values(ss) = (float)FS->ionSkinDepth[ss];
        saveToHDF5(group_scales, name, &fvec_values);
        #endif
        name.clear();

        delete group_scales;

        //Geometry of the mesh
        Group * group_geo = new Group( outputFile->createGroup( "/geometry" ) );

        name = "SPLIT_DIRECTION";
        int_value = params->mesh.SPLIT_DIRECTION;
        saveToHDF5(group_geo, name, &int_value);
        name.clear();

        name = "DX";
        cpp_type_value = params->mesh.DX;
        saveToHDF5(group_geo, name, &cpp_type_value);
        name.clear();

        name = "DY";
        cpp_type_value = params->mesh.DY;
        saveToHDF5(group_geo, name, &cpp_type_value);
        name.clear();

        name = "DZ";
        cpp_type_value = params->mesh.DZ;
        saveToHDF5(group_geo, name, &cpp_type_value);
        name.clear();

        name = "LX";
        cpp_type_value = params->mesh.LX;
        saveToHDF5(group_geo, name, &cpp_type_value);
        name.clear();

        name = "LY";
        cpp_type_value = params->mesh.LY;
        saveToHDF5(group_geo, name, &cpp_type_value);
        name.clear();

        name = "LZ";
        cpp_type_value = params->mesh.LZ;
        saveToHDF5(group_geo, name, &cpp_type_value);
        name.clear();

        name = "NX";
        int_value = params->mesh.NX_PER_MPI;
        saveToHDF5(group_geo, name, &int_value);
        name.clear();

        name = "NY";
        int_value = params->mesh.NY_PER_MPI;
        saveToHDF5(group_geo, name, &int_value);
        name.clear();

        name = "NZ";
        int_value = params->mesh.NZ_PER_MPI;
        saveToHDF5(group_geo, name, &int_value);
        name.clear();

        name = "NX_IN_SIM";
        int_value = params->mesh.NX_IN_SIM;
        saveToHDF5(group_geo, name, &int_value);
        name.clear();

        name = "NY_IN_SIM";
        int_value = params->mesh.NY_IN_SIM;
        saveToHDF5(group_geo, name, &int_value);
        name.clear();

        name = "NZ_IN_SIM";
        int_value = params->mesh.NZ_IN_SIM;
        saveToHDF5(group_geo, name, &int_value);
        name.clear();

        name = "X_MPI_CART_COORD";
        int_value = params->mpi.MPI_CART_COORDS_1D[0];
        saveToHDF5(group_geo, name, &int_value);
        name.clear();

        /*
        if (params->dimensionality == 1)
        {
          name = "X_MPI_CART_COORD";
          int_value = params->mpi.MPI_CART_COORDS_1D[0];
          saveToHDF5(group_geo, name, &int_value);
          name.clear();
        }
        else
        {
          name = "X_MPI_CART_COORD";
          int_value = params->mpi.MPI_CART_COORDS_2D[0];
          saveToHDF5(group_geo, name, &int_value);
          name.clear();

          name = "Y_MPI_CART_COORD";
          int_value = params->mpi.MPI_CART_COORDS_2D[1];
          saveToHDF5(group_geo, name, &int_value);
          name.clear();
        }
        */

        //Saving the x-axis coordinates
        name = "xAxis";

        #ifdef HDF5_DOUBLE
        vec_values = params->mesh.nodesX;
        saveToHDF5(group_geo, name, &vec_values);
        #elif defined HDF5_FLOAT
        fvec_values = conv_to<fvec>::from(params->mesh.nodesX);
        saveToHDF5(group_geo, name, &fvec_values);
        #endif

        name.clear();

        /*
        if (params->dimensionality == 2)
        {
          name = "yAxis";

          #ifdef HDF5_DOUBLE
          vec_values = params->mesh.nodesY;
          saveToHDF5(group_geo, name, &vec_values);
          #elif defined HDF5_FLOAT
          fvec_values = conv_to<fvec>::from(params->mesh.nodesY);
          saveToHDF5(group_geo, name, &fvec_values);
          #endif

          name.clear();
        }
        */

        delete group_geo;
        //Geometry of the mesh

        //Electron temperature
        name = "Te";
        cpp_type_value = params->f_IC.Te*F_KB/F_E;
        saveToHDF5(outputFile, name, &cpp_type_value);
        name.clear();

        //Ions
        Group * group_ions = new Group( outputFile->createGroup( "/ions" ) );

        name = "numberOfParticleSpecies";
        int_value = params->numberOfParticleSpecies;
        saveToHDF5(group_ions, name, &int_value);
        name.clear();

        name = "ne";
        cpp_type_value = (CPP_TYPE)params->f_IC.ne;
        saveToHDF5(group_ions, name, &cpp_type_value);
        name.clear();

        for(int ii=0;ii<params->numberOfParticleSpecies;ii++)
        {
          stringstream ionSpec;
          ionSpec << (ii+1);
          name = "/ions/spp_" + ionSpec.str();
          Group * group_ionSpecies = new Group( outputFile->createGroup( name ) );
          name.clear();

          name = "densityFraction";
          cpp_type_value = (CPP_TYPE)IONS->at(ii).p_IC.densityFraction;
          saveToHDF5(group_ionSpecies, name, &cpp_type_value);
          name.clear();

          name = "NCP";
          cpp_type_value = (CPP_TYPE)IONS->at(ii).NCP;
          saveToHDF5(group_ionSpecies, name, &cpp_type_value);
          name.clear();

          name = "NSP";
          cpp_type_value = (CPP_TYPE)IONS->at(ii).NSP;
          saveToHDF5(group_ionSpecies, name, &cpp_type_value);
          name.clear();

          name = "NSP_OUT";
          cpp_type_value = (CPP_TYPE)IONS->at(ii).nSupPartOutput;
          saveToHDF5(group_ionSpecies, name, &cpp_type_value);
          name.clear();

          name = "Tpar";
          cpp_type_value = (CPP_TYPE)IONS->at(ii).p_IC.Tpar*F_KB/F_E;
          saveToHDF5(group_ionSpecies, name, &cpp_type_value);
          name.clear();

          name = "Tper";
          cpp_type_value = (CPP_TYPE)IONS->at(ii).p_IC.Tper*F_KB/F_E;
          saveToHDF5(group_ionSpecies, name, &cpp_type_value);
          name.clear();

          name = "M";
          cpp_type_value = (CPP_TYPE)IONS->at(ii).M;
          saveToHDF5(group_ionSpecies, name, &cpp_type_value);
          name.clear();

          name = "Q";
          cpp_type_value = (CPP_TYPE)IONS->at(ii).Q;
          saveToHDF5(group_ionSpecies, name, &cpp_type_value);
          name.clear();

          name = "Z";
          cpp_type_value = (CPP_TYPE)IONS->at(ii).Z;
          saveToHDF5(group_ionSpecies, name, &cpp_type_value);
          name.clear();

          delete group_ionSpecies;
        }

        delete group_ions;
        //Ions

        //Electromagnetic fields
        name = "BX_0";
        vector_values = { (CPP_TYPE)params->em_IC.BX , (CPP_TYPE)params->em_IC.BX , (CPP_TYPE)params->em_IC.BX };
        saveToHDF5(outputFile, name, &vector_values);
        name.clear();

        delete outputFile;

      }//End of try block

      // catch failure caused by the H5File operations:
      // =============================================
      catch( FileIException error )
      {
            	error.printErrorStack();
      }

      // catch failure caused by the DataSet operations:
      // ===============================================
      catch( DataSetIException error )
      {
            	error.printErrorStack();
      }

      // catch failure caused by the DataSpace operations:
      // ================================================
      catch( DataSpaceIException error )
      {
            	error.printErrorStack();
      }
    } // MPI-0
}

void HDF_TYP::saveOutputs(const params_TYP * params, const vector<ionSpecies_TYP> * IONS, fields_TYP * fields, const CS_TYP * CS, const int it, double totalTime)
{

	try
	{
		stringstream iteration;
		stringstream dn;

		string name, path;
		path = params->PATH + "/HDF5/";

		dn << params->mpi.COMM_RANK;

		H5std_string FILE_NAME;

		// Save particles data
		if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
		{
			name = path + "PARTICLES_FILE_" + dn.str() + ".h5";
		}
		else if (params->mpi.COMM_COLOR == FIELDS_MPI_COLOR)
		{
			name = path + "FIELDS_FILE_" + dn.str() + ".h5";
		}

		FILE_NAME = name;

		name.clear();

		H5File * outputFile;

		if(it == 0)
		{
			outputFile = new H5File( FILE_NAME, H5F_ACC_TRUNC );// Create a new file using the default property lists.
		}
		else
		{
			outputFile = new H5File( FILE_NAME, H5F_ACC_RDWR );// Create a new file using the default property lists.
		}

		iteration << it;

		string group_iteration_name;
		group_iteration_name = "/" + iteration.str();
		Group * group_iteration = new Group( outputFile->createGroup( group_iteration_name ) );

		CPP_TYPE cpp_type_value;

		name = "time";
		cpp_type_value = (CPP_TYPE)totalTime;
		saveToHDF5(group_iteration, name, &cpp_type_value);
		name.clear();

		// Save particles data
		if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
		{
			saveIonsVariables(params, IONS, CS, group_iteration);
		}
		else if (params->mpi.COMM_COLOR == FIELDS_MPI_COLOR)
		{
			saveFieldsVariables(params, fields, CS, group_iteration);
		}

		delete group_iteration;

		delete outputFile;
	}//End of try block

	// catch failure caused by the H5File operations
    catch( FileIException error )
	{
		error.printErrorStack();
    }

	// catch failure caused by the DataSet operations
    catch( DataSetIException error )
	{
		error.printErrorStack();
    }

	// catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
	{
		error.printErrorStack();
    }
}


void HDF_TYP::saveIonsVariables(const params_TYP * params, const vector<ionSpecies_TYP> * IONS, const CS_TYP * CS, const Group * group_iteration)
{
	unsigned int iIndex(params->mesh.NX_PER_MPI*params->mpi.MPI_DOMAIN_NUMBER_CART+1);
	unsigned int fIndex(params->mesh.NX_PER_MPI*(params->mpi.MPI_DOMAIN_NUMBER_CART+1));

	try
	{
		string name;

		int int_value;
		CPP_TYPE cpp_type_value;
		std::vector<CPP_TYPE> vector_values;

		arma::ivec ivec_values;

		arma::vec vec_values;
		arma::fvec fvec_values;

		arma::mat mat_values;
		arma::fmat fmat_values;

		//Ions
		name = "ions";
		Group * group_ions = new Group( group_iteration->createGroup( "ions" ) );
		name.clear();

		//Iterations over the ion species.
		for(int ii=0; ii<IONS->size(); ii++)
		{
			stringstream ionSpec;
			ionSpec << (ii+1);
			name = "spp_" + ionSpec.str();
			Group * group_ionSpecies = new Group( group_ions->createGroup( name ) );
			name.clear();

			for(int ov=0; ov<params->outputs_variables.size(); ov++)
			{
				if(params->outputs_variables.at(ov) == "X_p")
				{
					//Saving the x-axis coordinates
					name = "X_p";
					#ifdef HDF5_DOUBLE
					vec_values = CS->length*IONS->at(ii).X_p.col(0);
					saveToHDF5(group_ionSpecies, name, &vec_values);
					#elif defined HDF5_FLOAT
					fvec_values = conv_to<fvec>::from(CS->length*IONS->at(ii).X_p.col(0));
					saveToHDF5(group_ionSpecies, name, &fvec_values);
					#endif
					name.clear();
				 }
				 /*
				 else if(params->outputs_variables.at(ov) == "pCount")
				 {
					name = "pCount";
					ivec_values = IONS->at(ii).pCount;
					saveToHDF5(group_ionSpecies, name, &ivec_values);
					name.clear();
				}
				else if(params->outputs_variables.at(ov) == "eCount")
				{
					//Saving ions species density
					name = "eCount";
					#ifdef HDF5_DOUBLE
					vec_values = (CS->velocity*CS->velocity*CS->mass)*IONS->at(ii).eCount;
					saveToHDF5(group_ionSpecies, name, &vec_values);
					#elif defined HDF5_FLOAT
					fvec_values = conv_to<fvec>::from((CS->velocity*CS->velocity*CS->mass)*IONS->at(ii).eCount);
					saveToHDF5(group_ionSpecies, name, &fvec_values);
					#endif
					name.clear();
                }
				*/
				else if(params->outputs_variables.at(ov) == "V_p")
				{
					name = "V_p";
					#ifdef HDF5_DOUBLE
					mat_values = CS->velocity*IONS->at(ii).V_p;

					saveToHDF5(group_ionSpecies, name, &mat_values);
					#elif defined HDF5_FLOAT
					fmat_values = conv_to<fmat>::from(CS->velocity*IONS->at(ii).V_p);
					saveToHDF5(group_ionSpecies, name, &fmat_values);
					#endif
					name.clear();
				}
				if(params->outputs_variables.at(ov) == "mu_p")
				{
					//Saving the x-axis coordinates
					name = "mu_p";
					#ifdef HDF5_DOUBLE
					vec_values = IONS->at(ii).mu_p*CS->magneticMoment;
					saveToHDF5(group_ionSpecies, name, &vec_values);
					#elif defined HDF5_FLOAT
					fvec_values = conv_to<fvec>::from(IONS->at(ii).mu_p*CS->magneticMoment);
					saveToHDF5(group_ionSpecies, name, &fvec_values);
					#endif
					name.clear();
				 }
				else if(params->outputs_variables.at(ov) == "EX_p")
				{
					name = "EX_p";
					#ifdef HDF5_DOUBLE
					vec_values = CS->eField*IONS->at(ii).EX_p;
					saveToHDF5(group_ionSpecies, name, &vec_values);
					#elif defined HDF5_FLOAT
					fvec_values = conv_to<fvec>::from( CS->eField*IONS->at(ii).EX_p);
					saveToHDF5(group_ionSpecies, name, &fvec_values);
					#endif
					name.clear();
				}
				else if(params->outputs_variables.at(ov) == "BX_p")
				{
					name = "BX_p";
					#ifdef HDF5_DOUBLE
					vec_values = CS->bField*IONS->at(ii).BX_p;
					saveToHDF5(group_ionSpecies, name, &vec_values);
					#elif defined HDF5_FLOAT
					fvec_values = conv_to<fvec>::from( CS->bField*IONS->at(ii).BX_p);
					saveToHDF5(group_ionSpecies, name, &fvec_values);
					#endif
					name.clear();
				}
				if(params->outputs_variables.at(ov) == "n_p")
				{
					//Saving the x-axis coordinates
					name = "n_p";
					#ifdef HDF5_DOUBLE
					vec_values = IONS->at(ii).n_p/CS->length;
					saveToHDF5(group_ionSpecies, name, &vec_values);
					#elif defined HDF5_FLOAT
					fvec_values = conv_to<fvec>::from(IONS->at(ii).n_p)/CS->length;
					saveToHDF5(group_ionSpecies, name, &fvec_values);
					#endif
					name.clear();
				 }
				 if(params->outputs_variables.at(ov) == "nv_p")
				 {
					 //Saving the "X" component of ion flux density at particle positions:
					 name = "nv_p";
					 #ifdef HDF5_DOUBLE
					 vec_values = IONS->at(ii).nv_p*CS->velocity/CS->length;
					 saveToHDF5(group_ionSpecies, name, &vec_values);
					 #elif defined HDF5_FLOAT
					 fvec_values = conv_to<fvec>::from(IONS->at(ii).nv_p)*CS->velocity/CS->length;
					 saveToHDF5(group_ionSpecies, name, &fvec_values);
					 #endif
					 name.clear();
				  }
				  if(params->outputs_variables.at(ov) == "Tpar_p")
				  {
					  //Saving the "X" drift velocity at the particle positions:
					  name = "Tpar_p";
					  #ifdef HDF5_DOUBLE
					  vec_values = IONS->at(ii).Tpar_p*CS->temperature*F_KB/F_E;
					  saveToHDF5(group_ionSpecies, name, &vec_values);
					  #elif defined HDF5_FLOAT
					  fvec_values = conv_to<fvec>::from(IONS->at(ii).Tpar_p)*CS->temperature*F_KB/F_E;
					  saveToHDF5(group_ionSpecies, name, &fvec_values);
					  #endif
					  name.clear();
				   }
				   if(params->outputs_variables.at(ov) == "Tper_p")
				   {
					   //Saving the "X" drift velocity at the particle positions:
					   name = "Tper_p";
					   #ifdef HDF5_DOUBLE
					   vec_values = IONS->at(ii).Tper_p*CS->temperature*F_KB/F_E;
					   saveToHDF5(group_ionSpecies, name, &vec_values);
					   #elif defined HDF5_FLOAT
					   fvec_values = conv_to<fvec>::from(IONS->at(ii).Tper_p)*CS->temperature*F_KB/F_E;
					   saveToHDF5(group_ionSpecies, name, &fvec_values);
					   #endif
					   name.clear();
					}
				else if(params->outputs_variables.at(ov) == "n_m")
				{
					if (params->mpi.IS_PARTICLES_ROOT)
					{
						//Saving ions species density
						name = "n_m";
						#ifdef HDF5_DOUBLE
						vec_values = IONS->at(ii).n_m.subvec(1,params->mesh.NX_IN_SIM)/CS->volume;
						saveToHDF5(group_ionSpecies, name, &vec_values);
						#elif defined HDF5_FLOAT
						fvec_values = conv_to<fvec>::from(IONS->at(ii).n_m.subvec(1,params->mesh.NX_IN_SIM)/CS->volume);
						saveToHDF5(group_ionSpecies, name, &fvec_values);
						#endif
						name.clear();
					}

				}
				else if(params->outputs_variables.at(ov) == "Tpar_m")
				{
					if (params->mpi.IS_PARTICLES_ROOT)
					{
						//Saving ions species density
						name = "Tpar_m";
						#ifdef HDF5_DOUBLE
						vec_values = IONS->at(ii).Tpar_m.subvec(1,params->mesh.NX_IN_SIM)*CS->temperature*F_KB/F_E;
						saveToHDF5(group_ionSpecies, name, &vec_values);
						#elif defined HDF5_FLOAT
						fvec_values = conv_to<fvec>::from(IONS->at(ii).Tpar_m.subvec(1,params->mesh.NX_IN_SIM)*CS->temperature*F_KB/F_E);
						saveToHDF5(group_ionSpecies, name, &fvec_values);
						#endif
						name.clear();
					}

				}
				else if(params->outputs_variables.at(ov) == "Tper_m")
				{
					if (params->mpi.IS_PARTICLES_ROOT)
					{
						//Saving ions species density
						name = "Tper_m";
						#ifdef HDF5_DOUBLE
						vec_values = IONS->at(ii).Tper_m.subvec(1,params->mesh.NX_IN_SIM)*CS->temperature*F_KB/F_E;
						saveToHDF5(group_ionSpecies, name, &vec_values);
						#elif defined HDF5_FLOAT
						fvec_values = conv_to<fvec>::from(IONS->at(ii).Tper_m.subvec(1,params->mesh.NX_IN_SIM)*CS->temperature*F_KB/F_E);
						saveToHDF5(group_ionSpecies, name, &fvec_values);
						#endif
						name.clear();
					}

				}
				else if(params->outputs_variables.at(ov) == "mn")
				{
					//Saving ions species density
					name = "mn";
					ivec_values = IONS->at(ii).mn;
					saveToHDF5(group_ionSpecies, name, &ivec_values);
					name.clear();

				}
				else if(params->outputs_variables.at(ov) == "u_m")
				{
					if (params->mpi.IS_PARTICLES_ROOT)
					{
						Group * group_bulkVelocity = new Group( group_ionSpecies->createGroup( "u_m" ) );

						//x-component species bulk velocity
						name = "x";
						#ifdef HDF5_DOUBLE
						vec_values = CS->velocity*IONS->at(ii).nv_m.subvec(1,params->mesh.NX_IN_SIM)/IONS->at(ii).n_m.subvec(1,params->mesh.NX_IN_SIM);
						saveToHDF5(group_ionSpecies, name, &vec_values);
						#elif defined HDF5_FLOAT
						fvec_values = conv_to<fvec>::from(CS->velocity*IONS->at(ii).nv_m.subvec(1,params->mesh.NX_IN_SIM)/IONS->at(ii).n_m.subvec(1,params->mesh.NX_IN_SIM));
						saveToHDF5(group_bulkVelocity, name, &fvec_values);
						#endif
						name.clear();

						/*
						//x-component species bulk velocity
						name = "y";
						#ifdef HDF5_DOUBLE
						vec_values = CS->velocity*IONS->at(ii).nv.Y.subvec(1,params->mesh.NX_IN_SIM)/IONS->at(ii).n.subvec(1,params->mesh.NX_IN_SIM);
						saveToHDF5(group_ionSpecies, name, &vec_values);
						#elif defined HDF5_FLOAT
						fvec_values = conv_to<fvec>::from(CS->velocity*IONS->at(ii).nv.Y.subvec(1,params->mesh.NX_IN_SIM)/IONS->at(ii).n.subvec(1,params->mesh.NX_IN_SIM));
						saveToHDF5(group_bulkVelocity, name, &fvec_values);
						#endif
						name.clear();

						//x-component species bulk velocity
						name = "z";
						#ifdef HDF5_DOUBLE
						vec_values = CS->velocity*IONS->at(ii).nv.Z.subvec(1,params->mesh.NX_IN_SIM)/IONS->at(ii).n.subvec(1,params->mesh.NX_IN_SIM);
						saveToHDF5(group_ionSpecies, name, &vec_values);
						#elif defined HDF5_FLOAT
						fvec_values = conv_to<fvec>::from(CS->velocity*IONS->at(ii).nv.Z.subvec(1,params->mesh.NX_IN_SIM)/IONS->at(ii).n.subvec(1,params->mesh.NX_IN_SIM));
						saveToHDF5(group_bulkVelocity, name, &fvec_values);
						#endif
						name.clear();
						*/

						delete group_bulkVelocity;
					}
				}
			}

			delete group_ionSpecies;
		}//Iterations over the ion species.

		delete group_ions;
	}//End of try block


    catch( FileIException error ){// catch failure caused by the H5File operations
		error.printErrorStack();
    }

    catch( DataSetIException error ){// catch failure caused by the DataSet operations
		error.printErrorStack();
    }

    catch( DataSpaceIException error ){// catch failure caused by the DataSpace operations
		error.printErrorStack();
    }
}


void HDF_TYP::saveFieldsVariables(const params_TYP * params, fields_TYP * fields, const CS_TYP * CS, const Group * group_iteration)
{
	unsigned int iIndex(params->mesh.NX_PER_MPI*params->mpi.MPI_DOMAIN_NUMBER_CART+1);
	unsigned int fIndex(params->mesh.NX_PER_MPI*(params->mpi.MPI_DOMAIN_NUMBER_CART+1));

	try{
		string name;

		int int_value;
		CPP_TYPE cpp_type_value;
		std::vector<CPP_TYPE> vector_values;

		arma::vec vec_values;
		arma::fvec fvec_values;

		arma::mat mat_values;
		arma::fmat fmat_values;

		Group * group_fields = new Group( group_iteration->createGroup( "fields" ) );//Electromagnetic fields

		for(int ov=0; ov<params->outputs_variables.size(); ov++)
		{
			if(params->outputs_variables.at(ov) == "EX_m")
			{
				Group * group_field = new Group( group_fields->createGroup( "EX_m" ) );//Electric fields

				//x-component of electric field
				name = "x";
				#ifdef HDF5_DOUBLE
				vec_values = CS->eField*fields->EX_m.subvec(iIndex,fIndex);
				saveToHDF5(group_ionSpecies, name, &vec_values);
				#elif defined HDF5_FLOAT
				fvec_values = conv_to<fvec>::from( CS->eField*fields->EX_m.subvec(iIndex,fIndex) );
				saveToHDF5(group_field, name, &fvec_values);
				#endif
				name.clear();

				/*
				//y-component of electric field
				name = "y";
				#ifdef HDF5_DOUBLE
				vec_values = CS->eField*fields->E.Y.subvec(iIndex,fIndex);
				saveToHDF5(group_ionSpecies, name, &vec_values);
				#elif defined HDF5_FLOAT
				fvec_values = conv_to<fvec>::from( CS->eField*fields->E.Y.subvec(iIndex,fIndex) );
				saveToHDF5(group_field, name, &fvec_values);
				#endif
				name.clear();

				//z-component of electric field
				name = "z";
				#ifdef HDF5_DOUBLE
				vec_values = CS->eField*fields->E.Z.subvec(iIndex,fIndex);
				saveToHDF5(group_ionSpecies, name, &vec_values);
				#elif defined HDF5_FLOAT
				fvec_values = conv_to<fvec>::from( CS->eField*fields->E.Z.subvec(iIndex,fIndex) );
				saveToHDF5(group_field, name, &fvec_values);
				#endif
				name.clear();
				*/

				delete group_field;
			}
			if(params->outputs_variables.at(ov) == "BX_m")
			{
				Group * group_field = new Group( group_fields->createGroup( "BX_m" ) );//Electric fields

				//x-component of magnetic field
				name = "x";
				#ifdef HDF5_DOUBLE
				vec_values = CS->bField*fields->BX_m.subvec(iIndex,fIndex);
				saveToHDF5(group_ionSpecies, name, &vec_values);
				#elif defined HDF5_FLOAT
				fvec_values = conv_to<fvec>::from( CS->bField*fields->BX_m.subvec(iIndex,fIndex) );
				saveToHDF5(group_field, name, &fvec_values);
				#endif
				name.clear();

				/*
				//y-component of magnetic field
				name = "y";
				#ifdef HDF5_DOUBLE
				vec_values = CS->bField*fields->B.Y.subvec(iIndex,fIndex);
				saveToHDF5(group_ionSpecies, name, &vec_values);
				#elif defined HDF5_FLOAT
				fvec_values = conv_to<fvec>::from( CS->bField*fields->B.Y.subvec(iIndex,fIndex) );
				saveToHDF5(group_field, name, &fvec_values);
				#endif
				name.clear();

				//z-component of magnetic field
				name = "z";
				#ifdef HDF5_DOUBLE
				vec_values = CS->bField*fields->B.Z.subvec(iIndex,fIndex);
				saveToHDF5(group_ionSpecies, name, &vec_values);
				#elif defined HDF5_FLOAT
				fvec_values = conv_to<fvec>::from( CS->bField*fields->B.Z.subvec(iIndex,fIndex) );
				saveToHDF5(group_field, name, &fvec_values);
				#endif
				name.clear();
				*/

				delete group_field;
			}
			if(params->outputs_variables.at(ov) == "dBX_m")
			{
				Group * group_field = new Group( group_fields->createGroup( "dBX_m" ) );//Electric fields

				name = "x";
				#ifdef HDF5_DOUBLE
				vec_values = fields->dBX_m.subvec(iIndex,fIndex)*CS->bField/CS->length;
				saveToHDF5(group_ionSpecies, name, &vec_values);
				#elif defined HDF5_FLOAT
				fvec_values = conv_to<fvec>::from( fields->dBX_m.subvec(iIndex,fIndex)*CS->bField/CS->length );
				saveToHDF5(group_field, name, &fvec_values);
				#endif
				name.clear();

				delete group_field;
			}
			if(params->outputs_variables.at(ov) == "ddBX_m")
			{
				Group * group_field = new Group( group_fields->createGroup( "ddBX_m" ) );//Electric fields

				name = "x";
				#ifdef HDF5_DOUBLE
				vec_values = fields->ddBX_m.subvec(iIndex,fIndex)*CS->bField/pow(CS->length,2);
				saveToHDF5(group_ionSpecies, name, &vec_values);
				#elif defined HDF5_FLOAT
				fvec_values = conv_to<fvec>::from( fields->ddBX_m.subvec(iIndex,fIndex)*CS->bField/pow(CS->length,2) );
				saveToHDF5(group_field, name, &fvec_values);
				#endif
				name.clear();

				delete group_field;
			}
		}

		delete group_fields;//Electromagnetic fields
	}


    catch( FileIException error ){// catch failure caused by the H5File operations
		error.printErrorStack();
    }

    catch( DataSetIException error ){// catch failure caused by the DataSet operations
		error.printErrorStack();
    }

    catch( DataSpaceIException error ){// catch failure caused by the DataSpace operations
		error.printErrorStack();
    }
}
