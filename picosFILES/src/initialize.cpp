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
map<string,string> init_TYP::readTextFile(string * inputFile)
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

// Constructor:
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
    params->errorCodes[-102] = "MPI's Cartesian topology could not be created";
    params->errorCodes[-103] = "Grid size violates assumptions of hybrid model for the plasma -- DX smaller than the electron skind depth can not be resolved";
    params->errorCodes[-106] = "Inconsistency in iniital ion's velocity distribution function";

    // Program information:
    // ===========================
    if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
        cout << "* PICOS++, a 1D-2V GC hybrid PIC code for Open plasma Systems           *" << endl;
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

}

// Populate params with data from input file:
// =============================================================================
void init_TYP::readInputFile(params_TYP * params)
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // ==================
    if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n";
        cout << "READING INPUT FILE ..." << endl;
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
	parametersStringMap = readTextFile(&name);

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
    params->mpi.MPIS_FIELDS  = stoi( parametersStringMap["mpisForFields"] );

    if(stoi( parametersStringMap["quietStart"] ) == 1)
    {
        params->quietStart = true;
    }
    else
    {
        params->quietStart = false;
    }

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
    //unsigned int NX      = (unsigned int)stoi( parametersStringMap["NX"] );
    params->dp           = stod( parametersStringMap["dp"] );
    params->geometry.r1  = stod( parametersStringMap["r1"] );
    params->geometry.r2  = stod( parametersStringMap["r2"] );
    params->geometry.LX_min  = stod( parametersStringMap["LX_min"] );
    params->geometry.LX_max  = stod( parametersStringMap["LX_max"] );

    // Electron initial conditions:
    // -------------------------------------------------------------------------
    params->f_IC.ne          = stod( parametersStringMap["IC_ne"] );
    params->f_IC.Te          = stod( parametersStringMap["IC_Te"] )*F_E/F_KB; // Te in eV in input file
    params->f_IC.Te_NX       = stoi( parametersStringMap["IC_Te_NX"] );
    params->f_IC.Te_fileName = parametersStringMap["IC_Te_fileName"];

    // RF parameters
    // -------------------------------------------------------------------------
    params->RF.Prf        = stod( parametersStringMap["RF_Prf"] );
    params->RF.n_harmonic = stoi( parametersStringMap["RF_n_harmonic"] );
    params->RF.freq       = stod( parametersStringMap["RF_freq"]);
    params->RF.x1         = stod( parametersStringMap["RF_x1"]  );
    params->RF.x2         = stod( parametersStringMap["RF_x2"]  );
    params->RF.t_ON       = stod( parametersStringMap["RF_t_ON"]  );
    params->RF.t_OFF      = stod( parametersStringMap["RF_t_OFF"]  );
    params->RF.kpar       = stod( parametersStringMap["RF_kpar"]);
    params->RF.kper       = stod( parametersStringMap["RF_kper"]);
    params->RF.handedness = stoi( parametersStringMap["RF_handedness"]);
    params->RF.Prf_NS     = stoi( parametersStringMap["RF_Prf_NS"] );
    params->RF.Prf_fileName = parametersStringMap["RF_Prf_fileName"];

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

    MPI_Barrier(MPI_COMM_WORLD);

    if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "READING INPUT FILE COMPLETED" << endl;
        cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n";
    }
}

// Read and populate ion parameters input file:
// =============================================================================
void init_TYP::readIonPropertiesFile(params_TYP * params, vector<ionSpecies_TYP> * IONS)
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
    parametersMap = readTextFile(&name);

    // Determine the total number of ION species:
    // =========================================
    int totalNumSpecies(params->numberOfParticleSpecies + params->numberOfTracerSpecies);

    // Loop over all ION species and extract data from "parametersMap":
    // ================================================================
    for(int ss=0;ss<totalNumSpecies;ss++)
    {
        string name;
        ionSpecies_TYP ions;
        int SPECIES;
        stringstream kk;

        kk << ss + 1;

        name = "SPECIES" + kk.str();
        SPECIES = stoi(parametersMap[name]);
        name.clear();

        if (SPECIES == 0 || SPECIES == 1)
        {
            // General:
            // =================================================================
            // Species type, 1: Guiding center, 0: Tracer
            ions.SPECIES = SPECIES;

            // Number of particles per cell:
            name = "NPC" + kk.str();
            ions.NPC = stoi(parametersMap[name]);
            name.clear();

            name = "pctSupPartOutput" + kk.str();
            ions.pctSupPartOutput = stod(parametersMap[name]);
            name.clear();

            // Charge state:
            name = "Z" + kk.str();
            ions.Z = stod(parametersMap[name]);
            name.clear();

            // AMU mass number:
            name = "M" + kk.str();
            ions.M = F_U*stod(parametersMap[name]);
            name.clear();

            // Initial condition:
            // =================================================================
            name = "IC_type_" + kk.str();
            ions.p_IC.IC_type = stoi(parametersMap[name]);
            name.clear();

            // Perpendicular temperature:
            name = "IC_Tper_" + kk.str();
            ions.p_IC.Tper = stod(parametersMap[name])*F_E/F_KB;
            name.clear();

            name = "IC_Tper_fileName_" + kk.str();
            ions.p_IC.Tper_fileName = parametersMap[name];
            name.clear();

            name = "IC_Tper_NX_" + kk.str();
            ions.p_IC.Tper_NX = stoi(parametersMap[name]);
            name.clear();

            // Parallel temperature:
            name = "IC_Tpar_" + kk.str();
            ions.p_IC.Tpar = stod(parametersMap[name])*F_E/F_KB; // Tpar in eV in input file
            name.clear();

            name = "IC_Tpar_fileName_" + kk.str();
            ions.p_IC.Tpar_fileName = parametersMap[name];
            name.clear();

            name = "IC_Tpar_NX_" + kk.str();
            ions.p_IC.Tpar_NX = stoi(parametersMap[name]);
            name.clear();

            // Density fraction relative to 1:
            name = "IC_densityFraction_" + kk.str();
            ions.p_IC.densityFraction = stod(parametersMap[name]);
            name.clear();

            name = "IC_densityFraction_fileName_" + kk.str();
            ions.p_IC.densityFraction_fileName = parametersMap[name];
            name.clear();

            name = "IC_densityFraction_NX_" + kk.str();
            ions.p_IC.densityFraction_NX = stoi(parametersMap[name]);
            name.clear();

            // Boundary conditions:
            // =================================================================
            name = "BC_type_" + kk.str();
            ions.p_BC.BC_type = stoi(parametersMap[name]);
            name.clear();

            name = "BC_T_" + kk.str();
            ions.p_BC.T = stod(parametersMap[name])*F_E/F_KB;
            name.clear();

            name = "BC_E_" + kk.str();
            ions.p_BC.E = stod(parametersMap[name])*F_E/F_KB;
            name.clear();

            name = "BC_eta_" + kk.str();
            ions.p_BC.eta = stod(parametersMap[name]);
            name.clear();

            name = "BC_G_" + kk.str();
            ions.p_BC.G = stod(parametersMap[name]);
            name.clear();

            name = "BC_mean_x_" + kk.str();
            ions.p_BC.mean_x = stod(parametersMap[name]);
            name.clear();

            name = "BC_sigma_x_" + kk.str();
            ions.p_BC.sigma_x = stod(parametersMap[name]);
            name.clear();

            name = "BC_G_fileName_" + kk.str();
            ions.p_BC.G_fileName = parametersMap[name];
            name.clear();

            name = "BC_G_NS_" + kk.str();
            ions.p_BC.G_NS = stoi(parametersMap[name]);
            name.clear();

            // Create new element on IONS vector:
            IONS->push_back(ions);
        }
        else
        {
            MPI_Barrier(MPI_COMM_WORLD);

            if(params->mpi.MPI_DOMAIN_NUMBER == 0)
            {
                cerr << "PRO++ ERROR: Enter a valid type of species -- options are 0 = tracers, 1 = guiding center" << endl;
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

// Read initial condition profiles from external files:
// =============================================================================
void init_TYP::readInitialConditionProfiles(params_TYP * params, electrons_TYP * electrons, vector<ionSpecies_TYP> * IONS)
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // ==================
    if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "* * * * * * * * * * * *  LOADING INITIAL CONDITION PROFILES  * * * * * * * * * * * * * * * * * *" << endl;
    }

    // Assemble PATH:
    // ==============
    int nn = params->PATH.length();
    std::string fileName  = params->PATH.substr(0,nn-13);

    // Assemble file paths:
    // ====================
    std::string fileName5 = fileName + "/inputFiles/" + params->em_IC.BX_fileName;
    std::string fileName6 = fileName + "/inputFiles/" + params->f_IC.Te_fileName;

    // Get electron temperature and magnetic field initial conditon data:
    // =================================================================
    params->em_IC.Bx_profile.load(fileName5);
    params->f_IC.Te_profile.load(fileName6);

    // Electron temperature "x" scale:
    // ===============================
    int Te_NX = params->f_IC.Te_NX;

    // Rescale profiles:
    // ================
    params->f_IC.Te_profile *= params->f_IC.Te;
    params->em_IC.Bx_profile *= params->em_IC.BX;

    // Get Ion initial condition profiles:
    // ==============================
    for(int ss=0;ss<IONS->size();ss++)
    {
        // Assemble  external filenames:
        // =============================
        std::string fileName2 = fileName + "/inputFiles/" + IONS->at(ss).p_IC.Tper_fileName;
        std::string fileName3 = fileName + "/inputFiles/" + IONS->at(ss).p_IC.Tpar_fileName;
        std::string fileName4 = fileName + "/inputFiles/" + IONS->at(ss).p_IC.densityFraction_fileName;

        // Load data from external file:
        // ============================
        IONS->at(ss).p_IC.Tper_profile.load(fileName2);
        IONS->at(ss).p_IC.Tpar_profile.load(fileName3);
        IONS->at(ss).p_IC.densityFraction_profile.load(fileName4);

        // Rescale the plasma Profiles:
        // ================================
        IONS->at(ss).p_IC.Tper_profile *= IONS->at(ss).p_IC.Tper;
        IONS->at(ss).p_IC.Tpar_profile *= IONS->at(ss).p_IC.Tpar;
        IONS->at(ss).p_IC.densityFraction_profile *= params->f_IC.ne;

        // Axial coordinate:
        // =================
        int Tper_NX = IONS->at(ss).p_IC.Tper_NX;
        IONS->at(ss).p_IC.x_profile = linspace(params->geometry.LX_min,params->geometry.LX_max,Tper_NX);
    }

    // Print to terminal:
    // ==================
    if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "* * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *" << endl;
        cout << "+ Succesfully read IC external files                        " << endl;
        cout << "* * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *" << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // ==================
    if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "* * * * * * * * * * * *  PLASMA INITIAL CONDITION PROFILES LOADED  * * * * * * * * * * * * * * * * * *" << endl;
    }

}

// Calculate IONS and params derived quantities:
// =============================================================================
void init_TYP::calculateDerivedQuantities(params_TYP * params, vector<ionSpecies_TYP> * IONS)
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // ==================
    if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << endl << "* * * * * * * * * * * CALCULATING DERIVED QUANTITIES * * * * * * * * * * * * * * * * *\n";
    }

    // Derived quantities for params:
    // ==============================
    params->mpi.MPIS_PARTICLES = params->mpi.NUMBER_MPI_DOMAINS - params->mpi.MPIS_FIELDS;
    params->geometry.A_0 = M_PI*(pow(params->geometry.r2,2) -pow(params->geometry.r1,2));
    params->geometry.LX = params->geometry.LX_max - params->geometry.LX_min;

    // Derived quantities for IONS:
    // ============================
    // NR total number of real particles at t = 0:
    double ds   = params->geometry.LX/params->em_IC.BX_NX;
    arma::vec A = params->geometry.A_0*(params->em_IC.BX/params->em_IC.Bx_profile);
    double ne   = params->CV.ne;
    double ne0  = params->f_IC.ne;
    double B    = params->CV.B;


    for(int ss=0; ss<IONS->size(); ss++)
    {
        // Ion species density:
        arma::vec n_ion = ones<vec>(params->em_IC.BX_NX);

        // Density fraction:
        double f = IONS->at(ss).p_IC.densityFraction;

        // Ion parameters:
        double M       = IONS->at(ss).M;
        double Z       = IONS->at(ss).Z;
        double Q       = F_E*Z;

        // Select ion density profile:
        if (params->quietStart)
        {
            n_ion *= ne0*f;
        }
        else
        {
            n_ion = IONS->at(ss).p_IC.densityFraction_profile;
        }

        // Number of real particles represented by particles of "ss" species:
        IONS->at(ss).NR = sum(n_ion%A)*ds;

        // Characteristic frequencies:
        IONS->at(ss).Q     = Q;
        IONS->at(ss).Wc    = Q*B/M;
        IONS->at(ss).Wp    = sqrt(f*ne*Q*Q/(F_EPSILON*M));

        // Characteristic thermal velocities:
        IONS->at(ss).VTper = sqrt(2.0*F_KB*params->CV.Tper/M);
        IONS->at(ss).VTpar = sqrt(2.0*F_KB*params->CV.Tpar/M);

        // Characteristic lengths and time scales:
        IONS->at(ss).LarmorRadius = IONS->at(ss).VTper/IONS->at(ss).Wc;
        IONS->at(ss).SkinDepth    = F_C/IONS->at(ss).Wp;
        IONS->at(ss).GyroPeriod   = 2.0*M_PI/IONS->at(ss).Wc;
    }

    // Characteristic length and time scales of simulation:
    // ===================================================
    // Assume that species 0 is the majority species
    params->ionLarmorRadius = IONS->at(0).LarmorRadius;
    params->ionSkinDepth    = IONS->at(0).SkinDepth;
    params->ionGyroPeriod   = IONS->at(0).GyroPeriod;

    // RF start and end time:
    // ======================
    params->RF.t_ON  *= params->ionGyroPeriod;
    params->RF.t_OFF *= params->ionGyroPeriod;

    // Estimate DX and NX:
    // ===================
    params->geometry.DX = params->dp*params->ionSkinDepth;
    params->geometry.NX = round(params->geometry.LX/params->geometry.DX);

    // Define final NX by making it a multiple of 2:
    // =============================================
    if (fmod(params->geometry.NX, 2.0) > 0.0)
    {
        params->geometry.NX -= 1;
    }

    // Define final value of "dp" and "DX":
    // ======================================
    params->dp = params->geometry.LX/(params->geometry.NX*params->ionSkinDepth);
    params->geometry.DX = params->dp*params->ionSkinDepth;

    // Define number of particles per MPI process:
    // ==========================================
    for(int ss=0; ss<IONS->size(); ss++)
    {
        // Number of computational particles for species "ss" per MPI process:
        IONS->at(ss).NSP = ceil( IONS->at(ss).NPC*(double)params->geometry.NX/(double)params->mpi.MPIS_PARTICLES );

        // Fraction of super particles to save in output:
        IONS->at(ss).nSupPartOutput = floor( (IONS->at(ss).pctSupPartOutput/100.0)*IONS->at(ss).NSP );

        // NCP conversion factor from super-particle to real particle:
        IONS->at(ss).NCP = IONS->at(ss).NR/(IONS->at(ss).NSP*params->mpi.MPIS_PARTICLES);
    }

    // Print main ion properties:
    // ===========================
    for(int ss=0; ss<IONS->size(); ss++)
    {
        if(params->mpi.MPI_DOMAIN_NUMBER == 0)
        {
            if (IONS->at(ss).SPECIES == 0)
            {
                cout << endl << "Species No "  << ss + 1 << " are tracers with the following parameters:" << endl;
            }
            else
            {
                cout << endl << "Species No "  << ss + 1 << " are guiding-center with the following parameters:" << endl;
            }

            // Stream to terminal:
            // ===================
            cout << "+ Number of real particles simulated at t = 0 s : " << IONS->at(ss).NR << endl;
            cout << "+ Number of real particles for every super-particle at t = 0 s : " << IONS->at(ss).NCP << endl;
            cout << "+ User-defined number of particles per MPI: " << IONS->at(ss).NSP << endl;
            cout << "+ Atomic number: " << IONS->at(ss).Z << endl;
            cout << "+ Mass: " << IONS->at(ss).M << " kg" << endl;
            cout << "+ Parallel temperature: " << params->CV.Tpar*F_KB/F_E << " eV" << endl;
            cout << "+ Perpendicular temperature: " << params->CV.Tper*F_KB/F_E << " eV" << endl;
            cout << "+ Cyclotron frequency: " << IONS->at(ss).Wc << " Hz" << endl;
            cout << "+ Plasma frequency: " << IONS->at(ss).Wp << " Hz" << endl;
            cout << "+ Parallel thermal velocity: " << IONS->at(ss).VTpar << " m/s" << endl;
            cout << "+ Perpendicular thermal velocity: " << IONS->at(ss).VTper << " m/s" << endl;
            cout << "+ Larmor radius: " << IONS->at(ss).LarmorRadius << " m" << endl;
        }
    }

    // Print to terminal:
    // ==================
    if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << endl << "* * * * * * * * * * * DERIVED QUANTITIES CALCULATION COMPLETE * * * * * * * * * * * * * * * * *\n";
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

// Calculate mesh geometry:
// =============================================================================
void init_TYP::calculateMeshParams(params_TYP * params)
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // ==================
    if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << endl << "* * * * * * * * * * * * LOADING/COMPUTING SIMULATION GRID * * * * * * * * * * * * * * * * * *\n";
    }

    params->mesh.LX         = params->geometry.LX;
    params->mesh.DX         = params->geometry.DX;
    params->mesh.NX_IN_SIM  = params->geometry.NX;
    params->mesh.NX_PER_MPI = (int)( (double)params->geometry.NX/(double)params->mpi.MPIS_FIELDS );

    // Set size of mesh and allocate memory:
    // =====================================
    params->mesh.nodesX.set_size(params->mesh.NX_IN_SIM);

    // Create mesh nodes: X domain
    // ============================
    for(int ii=0; ii<params->mesh.NX_IN_SIM; ii++)
    {
        // Mesh points at the center of each cell:
        params->mesh.nodesX(ii) = (double)ii*params->mesh.DX + (0.5*params->mesh.DX) + params->geometry.LX_min;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // ==================
    if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "+ Fraction of cell resolved: " << params->dp << endl;
        cout << "+ Number of mesh nodes along x-axis: " << params->mesh.NX_IN_SIM << endl;
    	cout << "+ Size of simulation domain along the x-axis: " << params->mesh.LX << " m" << endl;
    	cout << "* * * * * * * * * * * *  SIMULATION GRID LOADED/COMPUTED  * * * * * * * * * * * * * * * * * *" << endl;
    }
}

// Allocate memory to ION arrays:
// =============================================================================
void init_TYP::allocateMemoryIons(params_TYP * params, vector<ionSpecies_TYP> * IONS)
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Define total number of ions species:
    // ====================================
    int totalNumSpecies(params->numberOfParticleSpecies + params->numberOfTracerSpecies);

    // Print to terminal:
    // ==================
    if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << endl << "* * * * * * * * * * * * ALLOCATING MEMORY TO ION ARRAYS * * * * * * * * * * * * * * * * * *" << endl;
    }

    // Loop over all species:
    // ======================
    for (int ss=0; ss<totalNumSpecies; ss++)
    {
        allocateMeshDefinedIonArrays(params, &IONS->at(ss));

        if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
        {
            allocateParticleDefinedIonArrays(params, &IONS->at(ss));
        }

    }//Iteration over ion species

    // Print to terminal:
    // ==================
    if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "* * * * * * * * * * * * * MEMORY ALLOCATION COMPLETE P * * * * * * * * * * * * * * * * * * *" << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

// Allocate memory to mesh-defined ION arrays:
// =============================================================================
void init_TYP::allocateMeshDefinedIonArrays(const params_TYP * params, ionSpecies_TYP * IONS)
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

// Allocate memory to Particle-defined ION arrays:
// =============================================================================
void init_TYP::allocateParticleDefinedIonArrays(const params_TYP * params, ionSpecies_TYP * IONS)
{
    // Initialize particle-defined quantities:
    // ==================================
    IONS->X_p.zeros(IONS->NSP);
    IONS->V_p.zeros(IONS->NSP,2);
    IONS->mn.zeros(IONS->NSP);

    IONS->EX_p.zeros(IONS->NSP);
    IONS->BX_p.zeros(IONS->NSP);
    IONS->dBX_p.zeros(IONS->NSP);
    IONS->ddBX_p.zeros(IONS->NSP);

    IONS->wxc.zeros(IONS->NSP);
    IONS->wxl.zeros(IONS->NSP);
    IONS->wxr.zeros(IONS->NSP);

    IONS->n_p.zeros(IONS->NSP);
    IONS->nv_p.zeros(IONS->NSP);
    IONS->Tpar_p.zeros(IONS->NSP);
    IONS->Tper_p.zeros(IONS->NSP);

    IONS->Te_p.zeros(IONS->NSP);

    // Initialize particle defined flags:
    // ==================================
    IONS->f1.zeros(IONS->NSP);
    IONS->f2.zeros(IONS->NSP);
    IONS->f3.zeros(IONS->NSP);
    //IONS->f4.zeros(IONS->NSP);
    IONS->f5.zeros(IONS->NSP);

    // Initialize particle kinetic energy at boundaries:
    // ================================================
    IONS->dE1.zeros(IONS->NSP);
    IONS->dE2.zeros(IONS->NSP);
    IONS->dE3.zeros(IONS->NSP);
    //IONS->dE4.zeros(IONS->NSP);
    IONS->dE5.zeros(IONS->NSP);

    // Initialize resonance number:
    // ============================
    IONS->resNum.zeros(IONS->NSP);
    IONS->resNum_.zeros(IONS->NSP);

    // Rf terms:
    // ========
    IONS->udErf.zeros(IONS->NSP);
    IONS->doppler.zeros(IONS->NSP);
    IONS->udE3.zeros(IONS->NSP);

    // Initialize particle weight:
    // ===========================
    IONS->a_p.ones(IONS->NSP);

    // Initialize magnetic moment:
    // ===========================
    IONS->mu_p.zeros(IONS->NSP);
}

// Initialize ION particle position and velocity vector:
// =============================================================================
void init_TYP::initializeIons(const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS)
{
    MPI_Barrier(MPI_COMM_WORLD);

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
    for (int ss=0; ss<totalNumSpecies; ss++)
    {
        if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
        {
            initDist_TYP initDist(params);

            switch (IONS->at(ss).p_IC.IC_type)
            {
                case(1):
                {
                    if (params->quietStart)
                    {
                        initDist.uniform_maxwellianDistribution(params, &IONS->at(ss));
                    }
                    else
                    {
                        initDist.nonuniform_maxwellianDistribution(params, &IONS->at(ss));
                    }
                    break;
                }
                default:
                {
                }
            } // switch
        }

        // Print to the terminal:
        if(params->mpi.MPI_DOMAIN_NUMBER == 0)
        {
            cout << "ION SPECIES: " << (ss + 1) << endl;

            if (params->quietStart)
            {
                cout << "+ Using quiet start: YES" << endl;
            }
            else
            {
                cout << "+ Using quiet start: NO" << endl;
            }

            cout << "+ Super-particles used in simulation: " << IONS->at(ss).NSP*params->mpi.MPIS_PARTICLES << endl;
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

// Initialize electrons with profile data:
// =============================================================================
void init_TYP::initializeElectrons(const params_TYP * params, const CS_TYP * CS, electrons_TYP * electrons)
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // ==================
	if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << endl << "* * * * * * * * * * * * INITIALIZING ELECTRON FLUID * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
    }

    // Number of mesh points with ghost cells included:
    int NX(params->mesh.NX_IN_SIM + 2);

    // Allocate memory to the mesh-defined electron temperature:
    electrons->Te_m.zeros(NX);

    //Interpolate at mesh points:
    // ==========================
    // Query points:
    arma::vec xq = zeros(NX);
    arma::vec yq = zeros(NX);
    for(int ii=0; ii<NX; ii++)
    {
        xq(ii) = (double)ii*params->mesh.DX - (0.5*params->mesh.DX) + params->geometry.LX_min;
    }

    // Sample points:
    int Te_NX  = params->f_IC.Te_NX;

    // Spatial increment for external data:
    double dX = params->mesh.LX/((double)(Te_NX - 2));
    arma::vec xt = zeros(Te_NX);
    arma::vec yt = zeros(Te_NX);
    for(int ii=0; ii<Te_NX; ii++)
    {
        xt(ii) = (double)ii*dX - (0.5*dX) + params->geometry.LX_min;
    }

    // Te profile:
    // ===========
    arma::vec Te = params->f_IC.Te_profile;
    yt = Te;
    interp1(xt,yt,xq,yq);
    electrons->Te_m = yq;

    // Print to terminal:
    if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "Completed initializing the mesh-defined electron temperature" << endl;
        cout << "+ Reference electron temperature: " << scientific << params->f_IC.Te*F_KB/F_E << fixed << " [eV] " << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // ==================
	if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "* * * * * * * * * * * * ELECTRON FLUID INITIALIZED  * * * * * * * * * * * * * * * * * *" << endl;
    }

}

// Initialize fields with profile data:
// =============================================================================
void init_TYP::initializeFields(params_TYP * params, fields_TYP * fields)
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // ==================
	if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << endl << "* * * * * * * * * * * * INITIALIZING ELECTROMAGNETIC FIELDS * * * * * * * * * * * * * * * * * *" << endl;
    }

    // Number of mesh points with ghost cells included:
    int NX(params->mesh.NX_IN_SIM + 2);

    // Allocate memory to fields:
    fields->zeros(NX);

    // Select how to initialize electromagnetic fields:
    if (params->em_IC.uniformBfield)
    {
        fields->BX_m.fill(params->em_IC.BX);
    }
    else
    {
        //Interpolate at mesh points:
        // ==========================
        // Query points:
        arma::vec xq = zeros(NX);
        arma::vec yq = zeros(NX);
        for(int ii=0; ii<NX; ii++)
        {
            xq(ii) = (double)ii*params->mesh.DX - (0.5*params->mesh.DX) + params->geometry.LX_min;
        }

        // Sample points:
        int BX_NX  = params->em_IC.BX_NX;

        // Spatial increment for external data:
        double dX = params->mesh.LX/((double)(BX_NX - 2));
        arma::vec xt = zeros(BX_NX);
        arma::vec yt = zeros(BX_NX);
        for(int ii=0; ii<BX_NX; ii++)
        {
            xt(ii) = (double)ii*dX - (0.5*dX) + params->geometry.LX_min;
        }

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

    // Print to terminal:
    if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "Initializing electromagnetic fields within simulation" << endl;
        cout << "+ Reference magnetic field along x-axis: " << scientific << params->em_IC.BX << fixed << " T" << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // ==================
	if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "* * * * * * * * * * * * ELECTROMAGNETIC FIELDS INITIALIZED  * * * * * * * * * * * * * * * * * *" << endl;
    }
}
