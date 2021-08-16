#include "initialize.h"
#include "initDistribution.h"

// Function to split strings:
// =============================================================================
vector<string> init_TYP::split(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);

        if (pos == string::npos)
        {
            pos = str.length();
        }

        string token = str.substr(prev, pos-prev);

        if (!token.empty())
        {
            tokens.push_back(token);
        }

        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());

    return tokens;
}

// Function to read and load data from inputfile.input:
// =============================================================================
map<string,string> init_TYP::ReadAndloadInputFile(string * inputFile)
{
    // Create stream object:
    // =====================
    fstream reader;

    // Create map object:
    // ==================
    std::map<string,string> readMap;

    // Open input file using reader object:
    // ====================================
    reader.open(inputFile->data(),ifstream::in);

    // Handle error:
    // =============
    if (!reader)
    {
        MPI_Barrier(MPI_COMM_WORLD);

    	cerr << "PRO++ ERROR: The input file couldn't be opened." << endl;
    	MPI_Abort(MPI_COMM_WORLD, -101);
    }

    // Parse through file:
    // ===================
    string lineContent;
    vector<string> keyValuePair;
    while ( reader.good() )
    {
        // Read entire line:
        getline(reader,lineContent);

        // Search for comment symbol:
        size_t commentCharPos = lineContent.find("//",0);

        // Check for comment symbol:
        if (commentCharPos == 0 || lineContent.empty())
        {
            // Skip line
        }
        else
        {
            // Get value pair:
            keyValuePair = split(lineContent," ");

            // Update map:
            readMap[ keyValuePair[0] ] = keyValuePair[1];
        }
    }

    // Close stream object:
    // ===================
    reader.close();

    // Return map:
    // ==========
    return readMap;
}

// Constructor: Populates "params" based in inputfile.input
// =============================================================================
init_TYP::init_TYP(params_TYP * params, int argc, char* argv[])
{
    // Get RANK and SIZE of nodes within COMM_WORLD:
    // =============================================
    MPI_Comm_size(MPI_COMM_WORLD, &params->mpi.NUMBER_MPI_DOMAINS);
    MPI_Comm_rank(MPI_COMM_WORLD, &params->mpi.MPI_DOMAIN_NUMBER);

    // Error codes:
    // ============
    params->errorCodes[-100] = "Odd number of MPI processes";
    params->errorCodes[-101] = "Input file could not be opened";
    params->errorCodes[-102] = "MPI's Cartesian topology could not be created";
    params->errorCodes[-103] = "Grid size violates assumptions of hybrid model for the plasma -- DX smaller than the electron skind depth can not be resolved";
    params->errorCodes[-104] = "Loading external electromagnetic fields not implemented yet";
    params->errorCodes[-105] = "Restart not implemented yet";
    params->errorCodes[-106] = "Inconsistency in iniital ion's velocity distribution function";
    params->errorCodes[-107] = "Inconsistency in iniital ion's spatial distribution function";
    params->errorCodes[-108] = "Non-finite value in meshNode";
    params->errorCodes[-109] = "Number of nodes in either direction of simulation domain need to be a multiple of 2";
    params->errorCodes[-110] = "Non finite values in Ex";
    params->errorCodes[-111] = "Non finite values in Ey";
    params->errorCodes[-112] = "Non finite values in Ez";
    params->errorCodes[-113] = "Non finite values in Bx";
    params->errorCodes[-114] = "Non finite values in By";
    params->errorCodes[-115] = "Non finite values in Bz";

    // Copyright and Licence Info:
    // ===========================
    if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
        cout << "* PICOS++, a 1D-2V GC hybrid PIC code for open plasma systems           *" << endl;
        cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
        cout << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Arguments and paths to main function:
    // =====================================
    params->PATH = argv[2];
	params->argc = argc;
	params->argv = argv;

    // Check number of MPI domains:
    // ============================
	if( fmod( (double)params->mpi.NUMBER_MPI_DOMAINS, 2.0 ) > 0.0 )
    {
        MPI_Barrier(MPI_COMM_WORLD);

		if(params->mpi.MPI_DOMAIN_NUMBER == 0)
        {
			cerr << "PICOS++ ERROR: The number of MPI processes must be an even number." << endl;
		}

		MPI_Abort(MPI_COMM_WORLD,-100);
	}

    // Stream date when simulation is started:
    // =======================================
    if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        time_t current_time = std::time(NULL);
        cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << endl;
        cout << "STARTING " << params->argv[1] << " SIMULATION ON: " << std::ctime(&current_time) << endl;
        cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << endl;
    }

    // Get name of path to input file:
    // ===============================
	string name;
	if(params->argc > 3)
    {
		string argv(params->argv[3]);
		name = "inputFiles/input_file_" + argv + ".input";
		params->PATH += "/" + argv;
	}
    else
    {
		name = "inputFiles/input_file.input";
		params->PATH += "/";
	}

    // Read input file and assemble map:
    // ================================
	std::map<string,string> parametersStringMap;
	parametersStringMap = ReadAndloadInputFile(&name);

    // Create HDF5 folders if they don't exist:
    // ========================================
	if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
		string mkdir_outputs_dir = "mkdir " + params->PATH;
		const char * sys = mkdir_outputs_dir.c_str();
		int rsys = system(sys);

		string mkdir_outputs_dir_HDF5 = mkdir_outputs_dir + "/HDF5";
		sys = mkdir_outputs_dir_HDF5.c_str();
		rsys = system(sys);
	}

    // Populate "params" with data from input file:
    // ============================================

    // Assign input data from map to "params":
    // -------------------------------------------------------------------------
    params->mpi.MPIS_FIELDS  = stoi( parametersStringMap["mpisForFields"] );

    if(stoi( parametersStringMap["quietStart"] ) == 1)
    {
        params->quietStart = true;
    }
    else
    {
        params->quietStart = false;
    }

    params->numberOfRKIterations    = stoi( parametersStringMap["numberOfRKIterations"] );
    params->numberOfParticleSpecies = stoi( parametersStringMap["numberOfParticleSpecies"] );
    params->numberOfTracerSpecies   = stoi( parametersStringMap["numberOfTracerSpecies"] );
    params->advanceParticleMethod   = stoi( parametersStringMap["advanceParticleMethod"] );

    // Characteristic values:
    // -------------------------------------------------------------------------
    params->CV.ne   = stod( parametersStringMap["CV_ne"] );
    params->CV.Te   = stod( parametersStringMap["CV_Te"] )*F_E/F_KB;
    params->CV.B    = stod( parametersStringMap["CV_B"] );
    params->CV.Tpar = stod( parametersStringMap["CV_Tpar"] )*F_E/F_KB;
    params->CV.Tper = stod( parametersStringMap["CV_Tper"] )*F_E/F_KB;

    // Simulation time:
    // -------------------------------------------------------------------------
    params->DTc            = stod( parametersStringMap["DTc"] );
    params->simulationTime = std::stod( parametersStringMap["simulationTime"] );

    // Switches:
    // -------------------------------------------------------------------------
    params->SW.EfieldSolve   = stoi( parametersStringMap["SW_EfieldSolve"] );
    params->SW.HallTermSolve = stoi( parametersStringMap["SW_HallTermSolve"] );
    params->SW.BfieldSolve   = stoi( parametersStringMap["SW_BfieldSolve"] );
    params->SW.Collisions    = stoi( parametersStringMap["SW_Collisions"] );
    params->SW.RFheating     = stoi( parametersStringMap["SW_RFheating"] );
    params->SW.advancePos    = stoi( parametersStringMap["SW_advancePos"] );
    params->SW.linearSolve   = stoi( parametersStringMap["SW_linearSolve"] );

    // Magnetic field initial conditions:
    // -------------------------------------------------------------------------
    params->em_IC.uniformBfield = stoi( parametersStringMap["IC_uniformBfield"] );
    params->em_IC.BX            = stod( parametersStringMap["IC_BX"] );
    params->em_IC.BY            = stod( parametersStringMap["IC_BY"] );
    params->em_IC.BZ            = stod( parametersStringMap["IC_BZ"] );
    params->em_IC.BX_NX         = stoi( parametersStringMap["IC_BX_NX"] );
    params->em_IC.BX_fileName   = parametersStringMap["IC_BX_fileName"];

    // Geometry:
    // -------------------------------------------------------------------------
    unsigned int NX      = (unsigned int)stoi( parametersStringMap["NX"] );
    unsigned int NY      = (unsigned int)stoi( parametersStringMap["NY"] );
    unsigned int NZ      = (unsigned int)stoi( parametersStringMap["NZ"] );
    params->DrL          = stod( parametersStringMap["DrL"] );
    params->dp           = stod( parametersStringMap["dp"] );
    params->geometry.r1  = stod( parametersStringMap["r1"] );
    params->geometry.r2  = stod( parametersStringMap["r2"] );

    // Electron initial conditions:
    // -------------------------------------------------------------------------
    params->f_IC.ne      = stod( parametersStringMap["IC_ne"] );
    params->f_IC.Te      = stod( parametersStringMap["IC_Te"] )*F_E/F_KB; // Te in eV in input file

    // Output variables:
    // -------------------------------------------------------------------------
    params->outputCadence           = stod( parametersStringMap["outputCadence"] );
    string nonparsed_variables_list = parametersStringMap["outputs_variables"].substr(1, parametersStringMap["outputs_variables"].length() - 2);
    params->outputs_variables       = split(nonparsed_variables_list,",");

    // Data smoothing:
    // -------------------------------------------------------------------------
    params->smoothingParameter        = stod( parametersStringMap["smoothingParameter"] );
    params->filtersPerIterationFields = stoi( parametersStringMap["filtersPerIterationFields"] );
    params->filtersPerIterationIons   = stoi( parametersStringMap["filtersPerIterationIons"] );

    // Derived parameters:
    // ===================
    params->mpi.MPIS_PARTICLES = params->mpi.NUMBER_MPI_DOMAINS - params->mpi.MPIS_FIELDS;
    params->geometry.A_0 = M_PI*(pow(params->geometry.r2,2) -pow(params->geometry.r1,2));

    // Sanity check: if NX and/or NY is not a multiple of 2, the simulation aborts:
    // ============================================================================
    if (fmod(NX, 2.0) > 0.0)
    {
        MPI_Barrier(params->mpi.MPI_TOPO);
        MPI_Abort(params->mpi.MPI_TOPO,-109);
    }

    params->mesh.NX_IN_SIM = NX;
    params->mesh.NX_PER_MPI = (int)( (double)NX/(double)params->mpi.MPIS_FIELDS );
    params->mesh.SPLIT_DIRECTION = 0;

    params->mesh.NY_IN_SIM = 1;
    params->mesh.NY_PER_MPI = 1;

    params->mesh.NZ_IN_SIM = 1;
    params->mesh.NZ_PER_MPI = 1;

    params->mesh.NUM_CELLS_PER_MPI = params->mesh.NX_PER_MPI;
    params->mesh.NUM_CELLS_IN_SIM = params->mesh.NX_IN_SIM;

}

void init_TYP::loadIonParameters(params_TYP * params, vector<ionSpecies_TYP> * IONS)
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // ==================
    if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "* * * * * * * * * * * * LOADING ION PARAMETERS * * * * * * * * * * * * * * * * * *\n";
    	cout << "+ Number of ion species: " << params->numberOfParticleSpecies << endl;
    	cout << "+ Number of tracer species: " << params->numberOfTracerSpecies << endl;
    }

    // Assemble path to "ion_properties.ion":
    // ======================================
    string name;
    if(params->argc > 3)
    {
    	string argv(params->argv[3]);
    	name = "inputFiles/ions_properties_" + argv + ".ion";
    }
    else
    {
    	name = "inputFiles/ions_properties.ion";
    }

    // Read data from "ion_properties.ion" into "parametersMap":
    // =========================================================
    std::map<string,string> parametersMap;
    parametersMap = ReadAndloadInputFile(&name);

    // Determine the total number of ION species:
    // =========================================
    int totalNumSpecies(params->numberOfParticleSpecies + params->numberOfTracerSpecies);

    // Loop over all ION species and extract data from "parametersMap":
    // ================================================================
    for(int ii=0;ii<totalNumSpecies;ii++)
    {
        string name;
        ionSpecies_TYP ions;
        int SPECIES;
        stringstream ss;

        ss << ii + 1;

        name = "SPECIES" + ss.str();
        SPECIES = stoi(parametersMap[name]);
        name.clear();

        if (SPECIES == 0 || SPECIES == 1)
        {
            // General:
            // =================================================================
            // Species type, -1: Guiding center, 0: Tracer, 1: Full orbit
            ions.SPECIES = SPECIES;

            // Number of particles per cell:
            name = "NPC" + ss.str();
            ions.NPC = stoi(parametersMap[name]);
            name.clear();

            name = "pctSupPartOutput" + ss.str();
            ions.pctSupPartOutput = stod(parametersMap[name]);
            name.clear();

            // Charge state:
            name = "Z" + ss.str();
            ions.Z = stod(parametersMap[name]);
            name.clear();

            // AMU mass number:
            name = "M" + ss.str();
            ions.M = F_U*stod(parametersMap[name]);
            name.clear();

            // Initial condition:
            // =================================================================
            name = "IC_type_" + ss.str();
            ions.p_IC.IC_type = stoi(parametersMap[name]);
            name.clear();

            // Perpendicular temperature:
            name = "IC_Tper_" + ss.str();
            ions.p_IC.Tper = stod(parametersMap[name])*F_E/F_KB;
            name.clear();

            name = "IC_Tper_fileName_" + ss.str();
            ions.p_IC.Tper_fileName = parametersMap[name];
            name.clear();

            name = "IC_Tper_NX_" + ss.str();
            ions.p_IC.Tper_NX = stoi(parametersMap[name]);
            name.clear();

            // Parallel temperature:
            name = "IC_Tpar_" + ss.str();
            ions.p_IC.Tpar = stod(parametersMap[name])*F_E/F_KB; // Tpar in eV in input file
            name.clear();

            name = "IC_Tpar_fileName_" + ss.str();
            ions.p_IC.Tpar_fileName = parametersMap[name];
            name.clear();

            name = "IC_Tpar_NX_" + ss.str();
            ions.p_IC.Tpar_NX = stoi(parametersMap[name]);
            name.clear();

            // Density fraction relative to 1:
            name = "IC_densityFraction_" + ss.str();
            ions.p_IC.densityFraction = stod(parametersMap[name]);
            name.clear();

            name = "IC_densityFraction_fileName_" + ss.str();
            ions.p_IC.densityFraction_fileName = parametersMap[name];
            name.clear();

            name = "IC_densityFraction_NX_" + ss.str();
            ions.p_IC.densityFraction_NX = stoi(parametersMap[name]);
            name.clear();

            // Boundary conditions:
            // =================================================================
            name = "BC_type_" + ss.str();
            ions.p_BC.BC_type = stoi(parametersMap[name]);
            name.clear();

            name = "BC_T_" + ss.str();
            ions.p_BC.T = stod(parametersMap[name])*F_E/F_KB;
            name.clear();

            name = "BC_E_" + ss.str();
            ions.p_BC.E = stod(parametersMap[name])*F_E/F_KB;
            name.clear();

            name = "BC_eta_" + ss.str();
            ions.p_BC.eta = stod(parametersMap[name]);
            name.clear();

            name = "BC_G_" + ss.str();
            ions.p_BC.G = stod(parametersMap[name]);
            name.clear();

            name = "BC_mean_x_" + ss.str();
            ions.p_BC.mean_x = stod(parametersMap[name]);
            name.clear();

            name = "BC_sigma_x_" + ss.str();
            ions.p_BC.sigma_x = stod(parametersMap[name]);
            name.clear();

            name = "BC_G_fileName_" + ss.str();
            ions.p_BC.G_fileName = parametersMap[name];
            name.clear();

            name = "BC_G_NS_" + ss.str();
            ions.p_BC.G_NS = stoi(parametersMap[name]);
            name.clear();

            // Derived quantities:
            // =================================================================
            ions.Q     = F_E*ions.Z;
            ions.Wc    = ions.Q*params->CV.B/ions.M;
            ions.Wp    = sqrt( ions.p_IC.densityFraction*params->CV.ne*ions.Q*ions.Q/(F_EPSILON*ions.M) );//Check the definition of the plasma freq for each species!
            ions.VTper = sqrt(2.0*F_KB*params->CV.Tper/ions.M);
            ions.VTpar = sqrt(2.0*F_KB*params->CV.Tpar/ions.M);
            ions.LarmorRadius = ions.VTper/ions.Wc;

            //Definition of the initial total number of superparticles for each species
            ions.NSP = ceil( ions.NPC*(double)params->mesh.NUM_CELLS_IN_SIM/(double)params->mpi.MPIS_PARTICLES );

            // ******************
            // What is this for?
            ions.nSupPartOutput = floor( (ions.pctSupPartOutput/100.0)*ions.NSP );

            // Create new element on IONS vector:
            IONS->push_back(ions);

            // Print ion parameters to terminal:
            // =================================
            if(params->mpi.MPI_DOMAIN_NUMBER == 0)
            {
                if (ions.SPECIES == 0)
                {
                    cout << endl << "Species No "  << ii + 1 << " are tracers with the following parameters:" << endl;
                }
                else
                {
                    cout << endl << "Species No "  << ii + 1 << " are full-orbit particles with the following parameters:" << endl;
                }

                // Stream to terminal:
                // ===================
                cout << "+ User-defined number of particles per MPI: " << ions.NSP << endl;
                cout << "+ Atomic number: " << ions.Z << endl;
                cout << "+ Mass: " << ions.M << " kg" << endl;
                cout << "+ Parallel temperature: " << params->CV.Tpar*F_KB/F_E << " eV" << endl;
                cout << "+ Perpendicular temperature: " << params->CV.Tper*F_KB/F_E << " eV" << endl;
                cout << "+ Cyclotron frequency: " << ions.Wc << " Hz" << endl;
                cout << "+ Plasma frequency: " << ions.Wp << " Hz" << endl;
                cout << "+ Parallel thermal velocity: " << ions.VTpar << " m/s" << endl;
                cout << "+ Perpendicular thermal velocity: " << ions.VTper << " m/s" << endl;
                cout << "+ Larmor radius: " << ions.LarmorRadius << " m" << endl;
            }
        }
        else
        {
            MPI_Barrier(MPI_COMM_WORLD);

            if(params->mpi.MPI_DOMAIN_NUMBER == 0)
            {
                cerr << "PRO++ ERROR: Enter a valid type of species -- options are 0 = tracers, 1 = full orbit, -1 = guiding center" << endl;
            }
            MPI_Abort(MPI_COMM_WORLD,-106);
        }

    }//Iteration over ion species

    // Print to terminal:
    // ==================
    if(params->mpi.MPI_DOMAIN_NUMBER == 0)
        {
            cout << "* * * * * * * * * * * * ION PARAMETERS LOADED * * * * * * * * * * * * * * * * * *\n";
        }

    MPI_Barrier(MPI_COMM_WORLD);
}


void init_TYP::loadMeshGeometry(params_TYP * params, FS_TYP * FS)
{

    MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // ==================
    if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << endl << "* * * * * * * * * * * * LOADING/COMPUTING SIMULATION GRID * * * * * * * * * * * * * * * * * *\n";
    }

    // Select grid size: based on Larmour radius or ion skin depth:
    // ============================================================
    if( (params->DrL > 0.0) && (params->dp < 0.0) )
    {
    	params->mesh.DX = params->DrL*params->ionLarmorRadius;
    	params->mesh.DY = params->mesh.DX;
    	params->mesh.DZ = params->mesh.DX;

    	if(params->mpi.MPI_DOMAIN_NUMBER == 0)
            {
                cout << "Using LARMOR RADIUS to set up simulation grid." << endl;
            }
    }
    else if( (params->DrL < 0.0) && (params->dp > 0.0) )
    {
        params->mesh.DX = params->dp*params->ionSkinDepth;
        params->mesh.DY = params->mesh.DX;
        params->mesh.DZ = params->mesh.DX;

        if(params->mpi.MPI_DOMAIN_NUMBER == 0)
        {
            cout << "Using ION SKIN DEPTH to set up simulation grid." << endl;
        }
    }

    // Set size of mesh and allocate memory:
    // =====================================
    params->mesh.nodesX.set_size(params->mesh.NX_IN_SIM);
    params->mesh.nodesY.set_size(params->mesh.NY_IN_SIM);
    params->mesh.nodesZ.set_size(params->mesh.NZ_IN_SIM);

    // Create mesh nodes: X domain
    // ============================
    for(int ii=0; ii<params->mesh.NX_IN_SIM; ii++)
    {
        //params->mesh.nodesX(ii) = (double)ii*params->mesh.DX;
        // ****** does the above expresion reequire + 0.5*DX as in the fortran code?
        // ****** Notice that the above rnage goes from ii = 0 - (NX_IN_SIM - 1)

        // Mesh points at the center of each cell:
        params->mesh.nodesX(ii) = (double)ii*params->mesh.DX + 0.5*params->mesh.DX;
    }

    // Create mesh nodes: Y domain
    // ===========================
    for(int ii=0; ii<params->mesh.NY_IN_SIM; ii++)
    {
        params->mesh.nodesY(ii) = (double)ii*params->mesh.DY;
    }

    // Create mesh nodes: Z domain
    // ===========================
    for(int ii=0; ii<params->mesh.NZ_IN_SIM; ii++)
    {
        params->mesh.nodesZ(ii) = (double)ii*params->mesh.DZ;
    }

    // Define total length of each mesh:
    // =================================
    params->mesh.LX = params->mesh.DX*params->mesh.NX_IN_SIM;
    params->mesh.LY = params->mesh.DY*params->mesh.NY_IN_SIM;
    params->mesh.LZ = params->mesh.DZ*params->mesh.NZ_IN_SIM;

    // Print to terminal:
    // ==================
    if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "+ Number of mesh nodes along x-axis: " << params->mesh.NX_IN_SIM << endl;
        cout << "+ Number of mesh nodes along y-axis: " << params->mesh.NY_IN_SIM << endl;
        cout << "+ Number of mesh nodes along z-axis: " << params->mesh.NZ_IN_SIM << endl;

    	cout << "+ Size of simulation domain along the x-axis: " << params->mesh.LX << " m" << endl;
    	cout << "+ Size of simulation domain along the y-axis: " << params->mesh.LY << " m" << endl;
    	cout << "+ Size of simulation domain along the z-axis: " << params->mesh.LZ << " m" << endl;
    	cout << "* * * * * * * * * * * *  SIMULATION GRID LOADED/COMPUTED  * * * * * * * * * * * * * * * * * *" << endl;
    }
}


void init_TYP::loadPlasmaProfiles(params_TYP * params, vector<ionSpecies_TYP> * IONS)
{
    // Define number of mesh points with Ghost cells included:
    // ======================================================
    int NX(params->mesh.NX_IN_SIM + 2);

    // Define number of elements in external profiles:
    // ===============================================
    int nn = params->PATH.length();
    std::string inputFilePath = params->PATH.substr(0,nn-13);

    for(int ii=0;ii<IONS->size();ii++)
    {

        // Assemble  external filenames
        std::string fileName  = params->PATH.substr(0,nn-13);
        std::string fileName2 = fileName + "/inputFiles/" + IONS->at(ii).p_IC.Tper_fileName;
        std::string fileName3 = fileName + "/inputFiles/" + IONS->at(ii).p_IC.Tpar_fileName;
        std::string fileName4 = fileName + "/inputFiles/" + IONS->at(ii).p_IC.densityFraction_fileName;

        // Load data from external file:
        // ============================
        IONS->at(ii).p_IC.Tper_profile.load(fileName2);
        IONS->at(ii).p_IC.Tpar_profile.load(fileName3);
        IONS->at(ii).p_IC.densityFraction_profile.load(fileName4);

        // Rescale the plasma Profiles:
        // ================================
        IONS->at(ii).p_IC.Tper_profile *= IONS->at(ii).p_IC.Tper;
        IONS->at(ii).p_IC.Tpar_profile *= IONS->at(ii).p_IC.Tpar;
        IONS->at(ii).p_IC.densityFraction_profile *= params->f_IC.ne;

        // Print to terminal:
        // ==================
        if(params->mpi.MPI_DOMAIN_NUMBER == 0)
        {
            cout << "* * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *" << endl;
            cout << "+ Succesfully read external files for particle IC           " << endl;
            cout << "* * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *" << endl;
        }

    }

}


void init_TYP::initializeFields(params_TYP * params, fields_TYP * fields)
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // ==================
	if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << endl << "* * * * * * * * * * * * INITIALIZING ELECTROMAGNETIC FIELDS * * * * * * * * * * * * * * * * * *" << endl;
    }

    // Initialize fields:
    // ================================
    initializeFieldsSizeAndValue(params,fields);

    // Print to terminal:
    if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "Initializing electromagnetic fields within simulation" << endl;
        cout << "+ Magnetic field along x-axis: " << scientific << params->em_IC.BX << fixed << " T" << endl;
        cout << "+ Magnetic field along y-axis: " << scientific << params->em_IC.BY << fixed << " T" << endl;
        cout << "+ Magnetic field along z-axis: " << scientific << params->em_IC.BZ << fixed << " T" << endl;
    }

    // Print to terminal:
    // ==================
	if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "* * * * * * * * * * * * ELECTROMAGNETIC FIELDS INITIALIZED  * * * * * * * * * * * * * * * * * *" << endl;
    }
}

void init_TYP::initializeFieldsSizeAndValue(params_TYP * params, fields_TYP * fields)
{

    // Number of mesh points with ghost cells included:
    // ================================================
    int NX(params->mesh.NX_IN_SIM + 2);

    // Allocate memory to fields:
    // ==========================
    fields->zeros(NX);

    // Select how to initialize electromagnetic fields:
    // ================================================
    if (params->em_IC.uniformBfield)
    {
        fields->BX_m.fill(params->em_IC.BX);
    }
    else
    {
        // Read filename from input file:
        // =============================
        int nn = params->PATH.length();
        std::string fileName = params->PATH.substr(0,nn-13);
        fileName = fileName + "/inputFiles/" + params->em_IC.BX_fileName;

        // Load data from external file:
        // ============================
        params->em_IC.Bx_profile.load(fileName);

        // Scale the normalize profile to Tesla:
        // ====================================
        params->em_IC.Bx_profile *= params->em_IC.BX;

        //Interpolate at mesh points:
        // ==========================
        // Query points:
        arma::vec xq = linspace(0,params->mesh.LX,NX);
        arma::vec yq(xq.size());

        // Sample points:
        int BX_NX  = params->em_IC.BX_NX;

        // Spatial increment for external data:
        double dX = params->mesh.LX/((double)BX_NX);

        arma::vec xt = linspace(0,params->mesh.LX,BX_NX); // x-vector from the table
        arma::vec yt(xt.size());

        // BX profile:
        // ===========
        arma::vec BX = params->em_IC.Bx_profile;
        yt = BX;
        interp1(xt,yt,xq,yq);
        fields->BX_m = yq;

        // dBX profile:
        // ===========
        arma::vec dBX(BX_NX,1);
        dBX.subvec(1,BX_NX-2) = (BX.subvec(2,BX_NX-1) - BX.subvec(0,BX_NX-3))/(2*dX);
        dBX(0)       = dBX(1);
        dBX(BX_NX-1) = dBX(BX_NX-2);
        yt = dBX;
        interp1(xt,yt,xq,yq);
        fields->dBX_m = yq;

        // ddBX profile:
        // ============
        arma::vec ddBX(BX_NX,1);
        ddBX.subvec(0,BX_NX-3) = diff(params->em_IC.Bx_profile,2)/(dX*dX);
        ddBX(BX_NX-2) = ddBX(BX_NX-3);
        ddBX(BX_NX-1) = ddBX(BX_NX-2);

        yt = ddBX;
        interp1(xt,yt,xq,yq);
        fields->ddBX_m = yq;
    }
}


void init_TYP::setupIonsInitialCondition(const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    // Define total number of ions species:
    // ====================================
    int totalNumSpecies(params->numberOfParticleSpecies + params->numberOfTracerSpecies);

    // Print to terminal:
    // ==================
    if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << endl << "* * * * * * * * * * * * SETTING UP IONS INITIAL CONDITION * * * * * * * * * * * * * * * * * *" << endl;
    }

    // Loop over all species:
    // ======================
    for (int ii=0; ii<totalNumSpecies; ii++)
    {
        if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
        {
            // Create start object:
            initDist_TYP initDist(params);

            switch (IONS->at(ii).p_IC.IC_type)
            {
                case(1):
                {
                    if (params->quietStart)
                    {
                        initDist.uniform_maxwellianDistribution(params, &IONS->at(ii));
                    }
                    else
                    {
                        initDist.nonuniform_maxwellianDistribution(params, &IONS->at(ii));
                    }
                    break;
                }
                default:
                {
                }
            } // switch

            initializeParticlesArrays(params, fields, &IONS->at(ii));

            initializeBulkVariablesArrays(params, &IONS->at(ii));
        }
        else if (params->mpi.COMM_COLOR == FIELDS_MPI_COLOR)
        {
            initializeBulkVariablesArrays(params, &IONS->at(ii));
        }

        // Broadcast NSP (Number of super particles per process) and nSupPartPutput from ROOTS to COMM_WORLD:
        MPI_Bcast(&IONS->at(ii).NSP, 1, MPI_DOUBLE, params->mpi.PARTICLES_ROOT_WORLD_RANK, MPI_COMM_WORLD);
        MPI_Bcast(&IONS->at(ii).nSupPartOutput, 1, MPI_DOUBLE, params->mpi.PARTICLES_ROOT_WORLD_RANK, MPI_COMM_WORLD);

        // NR total number of real particles:
        double ds    = params->mesh.LX/params->em_IC.BX_NX;
        arma::vec ni = ones<vec>(params->em_IC.BX_NX);
        arma::vec A  = ones<vec>(params->em_IC.BX_NX);
        if (params->quietStart)
        {
            ni *= params->f_IC.ne*IONS->at(ii).p_IC.densityFraction;
            A  = A*params->geometry.A_0;
        }
        else
        {
            ni = IONS->at(ii).p_IC.densityFraction_profile;
            A  = params->geometry.A_0*(params->em_IC.BX/params->em_IC.Bx_profile);
        }
        double NR    = sum(ni%A)*ds;

        // NCP conversion factor from super-particle to real particle:
        double NSP  = IONS->at(ii).NSP*params->mpi.MPIS_PARTICLES;
        IONS->at(ii).NCP = NR/NSP;

        if(params->mpi.MPI_DOMAIN_NUMBER == 0)
        {
            cout << "NR = " << NR << endl;
            cout << "NCP = " << IONS->at(ii).NCP << endl;
        }

        // Print to the terminal:
        if(params->mpi.MPI_DOMAIN_NUMBER == 0)
        {
            cout << "iON SPECIES: " << (ii + 1) << endl;

            if (params->quietStart)
            {
                cout << "+ Using quiet start: YES" << endl;
            }
            else
            {
                cout << "+ Using quiet start: NO" << endl;
            }

            cout << "+ Super-particles used in simulation: " << IONS->at(ii).NSP*params->mpi.MPIS_PARTICLES << endl;
        }

    }//Iteration over ion species

    // Print to terminal:
    // ==================
    if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "* * * * * * * * * * * * * IONS INITIAL CONDITION SET UP * * * * * * * * * * * * * * * * * * *" << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void init_TYP::initializeParticlesArrays(const params_TYP * params, fields_TYP * fields, ionSpecies_TYP * IONS)
{
    // Set size and value to zero of arrays for ions' variables:
    // ========================================================
    IONS->mn.zeros(IONS->NSP);

    IONS->EX_p.zeros(IONS->NSP);
    IONS->BX_p.zeros(IONS->NSP);
    IONS->dBX_p.zeros(IONS->NSP);
    IONS->ddBX_p.zeros(IONS->NSP);

    IONS->wxc.zeros(IONS->NSP);
    IONS->wxl.zeros(IONS->NSP);
    IONS->wxr.zeros(IONS->NSP);

    IONS->wxc_.zeros(IONS->NSP);
    IONS->wxl_.zeros(IONS->NSP);
    IONS->wxr_.zeros(IONS->NSP);

    // Initialize particle-defined quantities:
    // ==================================
    IONS->n_p.zeros(IONS->NSP);
    IONS->nv_p.zeros(IONS->NSP);
    IONS->Tpar_p.zeros(IONS->NSP);
    IONS->Tper_p.zeros(IONS->NSP);

    // Initialize particle defined flags:
    // ==================================
    IONS->f1.zeros(IONS->NSP);
    IONS->f2.zeros(IONS->NSP);
    IONS->f3.zeros(IONS->NSP);

    // Initialize particle kinetic energy at boundaries:
    // ================================================
    IONS->dE1.zeros(IONS->NSP);
    IONS->dE2.zeros(IONS->NSP);
    IONS->dE3.zeros(IONS->NSP);

    // Initialize particle weight:
    // ===========================
    IONS->a_p.ones(IONS->NSP);

    // Initialize magnetic moment:
    // ===========================
    IONS->mu_p.zeros(IONS->NSP);

    // Assign cell:
    // ===========
    // Populates wxc,wxl and wxr only
    //PIC_TYP PIC;
    //PIC.assignCell(params, fields, IONS);
}

void init_TYP::initializeBulkVariablesArrays(const params_TYP * params, ionSpecies_TYP * IONS)
{
    // Initialize mesh-defined quantities:
    // ==================================
    // Ion density:
    IONS->n_m.zeros(params->mesh.NX_IN_SIM + 2);
    IONS->n_m_.zeros(params->mesh.NX_IN_SIM + 2);
    IONS->n_m__.zeros(params->mesh.NX_IN_SIM + 2);
    IONS->n_m___.zeros(params->mesh.NX_IN_SIM + 2);

    // Ion flux density:
    IONS->nv_m.zeros(params->mesh.NX_IN_SIM + 2);
    IONS->nv_m_.zeros(params->mesh.NX_IN_SIM + 2);
    IONS->nv_m__.zeros(params->mesh.NX_IN_SIM + 2);

    // Pressure tensors:
    IONS->P11_m.zeros(params->mesh.NX_IN_SIM + 2);
    IONS->P22_m.zeros(params->mesh.NX_IN_SIM + 2);

    // Derived quantities:
    IONS->Tpar_m.zeros(params->mesh.NX_IN_SIM + 2);
    IONS->Tper_m.zeros(params->mesh.NX_IN_SIM + 2);
}
