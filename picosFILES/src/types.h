#ifndef H_TYPES
#define H_TYPES

#include <typeinfo>
#include <iostream>
#include <vector>
#include <armadillo>
#include <string>
#include <map>

#include<omp.h>
#include "mpi.h"

using namespace std;

//  Define macros:
// =============================================================================
#define FIELDS_MPI_COLOR 0
#define PARTICLES_MPI_COLOR 1
#define FIELDS_TAG 100
#define PARTICLES_TAG 200

#define float_zero 1E-7
#define double_zero 1E-15

// Physical constants
// =============================================================================
#define PRO_ZERO 1.0E-15	// Definition of zero in PROMETHEUS
#define F_E 1.602176E-19	// Electron charge in C (absolute value)
#define F_ME 9.109382E-31	// Electron mass in kg
#define F_MP 1.672621E-27	// Proton mass in kg
#define F_U 1.660538E-27	// Atomic mass unit in kg
#define F_KB 1.380650E-23	// Boltzmann constant in Joules/Kelvin
#define F_EPSILON 8.854E-12 // Vacuum permittivity in C^2/(N*m^2)
#define F_C 299792458.0 	// Light speed in m/s
#define F_MU (4*M_PI)*1E-7 	// Vacuum permeability in N/A^2
extern double F_EPSILON_DS; // Dimensionless vacuum permittivity
extern double F_E_DS; 		// Dimensionless electron charge
extern double F_ME_DS; 		// Dimensionless electron mass
extern double F_MU_DS; 		// Dimensionless vacuum permeability
extern double F_C_DS; 		// Dimensionless speed of light

// Declare vector field types:
// =============================================================================
class vfield_vec_TYP
{

public:
	arma::vec X;
	arma::vec Y;
	arma::vec Z;

	vfield_vec_TYP(){};

	vfield_vec_TYP(unsigned int N);

	~vfield_vec_TYP(){};

  /*
	vfield_vec_TYP operator + (vfield_vec_TYP R);
	vfield_vec_TYP operator += (vfield_vec_TYP R);
	vfield_vec_TYP operator - (vfield_vec_TYP R);
	vfield_vec_TYP operator -= (vfield_vec_TYP R);
	vfield_vec_TYP operator * (double s);
	vfield_vec_TYP operator *= (double s);
  friend vfield_vec_TYP operator * (double s, vfield_vec_TYP R );
	vfield_vec_TYP operator / (double s);
	vfield_vec_TYP operator / (vfield_vec_TYP R);
	vfield_vec_TYP operator /= (double s);
	vfield_vec_TYP operator /= (vfield_vec_TYP R);
  */

  /*
	void fill(double value);
	void ones(unsigned int N);
	void zeros();
	void zeros(unsigned int N);
  */

};

// Declare matrix field types:
// =============================================================================
class vfield_mat_TYP
{

public:
	arma::mat X;
	arma::mat Y;
	arma::mat Z;

	vfield_mat_TYP(){};

  /*
	vfield_mat_TYP(unsigned int N, unsigned int M);
  */

	~vfield_mat_TYP(){};

  /*
	vfield_mat_TYP operator + (vfield_mat_TYP R);
	vfield_mat_TYP operator += (vfield_mat_TYP R);
	vfield_mat_TYP operator - (vfield_mat_TYP R);
	vfield_mat_TYP operator -= (vfield_mat_TYP R);
	vfield_mat_TYP operator * (double s);
	vfield_mat_TYP operator *= (double s);
  friend vfield_mat_TYP operator * (double s, vfield_mat_TYP R );
	vfield_mat_TYP operator / (double s);
	vfield_mat_TYP operator / (vfield_mat_TYP R);
	vfield_mat_TYP operator /= (double s);
	vfield_mat_TYP operator /= (vfield_mat_TYP R);
  */

  /*
	void fill(double value);
	void ones(unsigned int N, unsigned int M);
	void zeros();
	void zeros(unsigned int N, unsigned int M);
  */

};

//  Structure to store each ion species initial condition parameters:
// =============================================================================
struct p_IC_TYP
{
	int IC_type;             	   		// 1: Uniform profiles, 2: profiles from external files

	// Reference values for profiles:
	// ==============================
	double Tper;
	double Tpar;
	double densityFraction;

	// Name of external files:
	// =======================
	string Tper_fileName; 				// File containing normalized spatial profile of Tper
	string Tpar_fileName; 				// File containing normalized spatial profile of Tpar
	string densityFraction_fileName;	// File containing normalized spatial profile

	// Number of elements of profiles:
	// ===============================
	int Tper_NX;
	int Tpar_NX;
	int densityFraction_NX;

	// Variables to store profiles from external files:
	// ================================================
	arma::vec Tper_profile;
	arma::vec Tpar_profile;
	arma::vec densityFraction_profile;
};

//  Structure to store each ion species particle boundary condition parameters:
// =============================================================================
struct p_BC_TYP
{
	int BC_type;					// 1: Warm plasma source, 2: NBI, 3: periodic, 4: reflecting

	// Particle source temperature:
	// ===========================
	double T;

	// Particle source beam energy:
	// ============================
	double E;

	// Particle source beam pitch angle:
	// ================================
	double eta;

	// Particle source rate:
	// =====================
	double G;						// Fueling rate of source in particles/second
	double sigma_x;					// Spatial spread of beam
	double mean_x;					// Spatial location of beam injection

	// Name of external file:
	// ======================
	string G_fileName;

	// Number of elements in external file:
	// ====================================
	int G_NS;

	// Variable to store profile from external file:
	// =============================================
	arma::vec G_profile;

	// Computational particle accumulators:
	// =============================================
	double S1;
	double S2;
	double GSUM;

	// Variable to store new particle weight based on fueling:
	// ======================================================
	double a_new;

	// Constructor:
	// ===========
	p_BC_TYP()
	{
		BC_type   = 0;
		E 	  	  = 0;
		T     	  = 0;
		eta   	  = 0;
		G         = 0;
		sigma_x   = 0;
		mean_x    = 0;
		G_NS      = 0;
		GSUM      = 0;
		S1        = 0;
		S2        = 0;
		a_new 	  = 0;
	};

};

//  Define ION VARIABLES AND PARAMETERS DERIVED TYPES:
// =============================================================================
class ionSpecies_TYP : public vfield_vec_TYP
{

public:

	int SPECIES;
	int IC; 					// Initial condition IC=1 (Maxwellian), IC=2 (ring-like)

	double NSP; 				// Initial number of superparticles for the given ion species.
	double NCP; 				// Number of charged particles per superparticle.
	double NPC; 				// Number of superparticles per cell. When its value is zero, the particles are loaded from external files.
	double Q; 					// Charge.
	double Z; 					// Atomic number.
	double M; 					// Mass

	// variables for controlling super-particles' outputs
	double pctSupPartOutput; 	//
	unsigned int nSupPartOutput;//

	//double densityFraction;

	double go;					// Initial relativistic gamma
	double LarmorRadius;		// Larmor radius.
	double VTper;				// Thermal velocity.
	double VTpar;				// Thermal velocity.
	double Wc;					// Average cyclotron frequency.
	double Wp;					// Plasma frequency.
	double avg_mu; 				// Average magnetic moment

	arma::mat X; 				// Ions position, the dimension should be (NSP,3), where NP is the number of particles of the ion species.
	arma::mat V; 				// Ions' velocity, the dimension should be (NSP,3), where NP is the number of particles of the ion species.
	arma::mat P; 				// Ions' momentum, the dimension should be (NSP,3), where NP is the number of particles of the ion species.
	arma::vec g; 				// Ions' relativistic gamma factor.
	arma::ivec mn; 			// Ions' position in terms of the index of mesh node


  // ************ Consider arma::vec E, arma::vec B, arma::vec dB , arma::vec ddB ***********
	arma::mat E;				// Electric field seen by particles when advancing particles velocity
	arma::mat B;				// Magnetic field seen by particles when advancing particles velocity
// ***********************


  // ***************** Consider deleting ************
	// Guiding-center variables
	//arma::vec mu; 				// Ions' magnetic moment.
	//arma::vec Ppar; 			// Parallel momentum used in guiding-center orbits
  // ***************** Consider deleting ************

	//These weights are used in the charge extrapolation and the force interpolation
	arma::vec wxl;				// Particles' weights w.r.t. the vertices of the grid cells
	arma::vec wxc;				// Particles' weights w.r.t. the vertices of the grid cells
	arma::vec wxr;				// Particles' weights w.r.t. the vertices of the grid cells

	arma::vec wxl_;				// Particles' weights w.r.t. the vertices of the grid cells
	arma::vec wxc_;				// Particles' weights w.r.t. the vertices of the grid cells
	arma::vec wxr_;				// Particles' weights w.r.t. the vertices of the grid cells

    // Mesh-defined ion moments:
	arma::vec n_m; 				// Ion density at time level "l + 1"
	arma::vec n_m_; 				// Ion density at time level "l"
	arma::vec n_m__; 				// Ion density at time level "l - 1"
	arma::vec n_m___; 			// Ion density at time level "l - 2"
	vfield_vec_TYP nv_m; 				// Ion bulk velocity at time level "l + 1/2"
	vfield_vec_TYP nv_m_; 			// Ion bulk velocity at time level "l - 1/2"
	vfield_vec_TYP nv_m__; 			// Ion bulk velocity at time level "l - 3/2"
  arma::vec P11_m;				// Ion pressure tensor, component 1,1
  arma::vec P22_m;				// Ion pressure tensor, component 2,2
  arma::vec Tpar_m;			// Ion parallel temperature
  arma::vec Tper_m;			// Ion perpendicular temperature

	// Particle-defined ion moments:
	arma::vec n_p;
	arma::vec nv_p;
	arma::vec Tpar_p;
	arma::vec Tper_p;

	// Particle defined flags:
	arma::ivec f1;             	// Flag for left boundary
	arma::ivec f2;              // Flag for Right boundary
	arma::ivec f3;              // Flag for RF operator

	// Particle kinetic energy at boundaries:
	arma::ivec dE1;              // left boundary
	arma::ivec dE2;              // Right boundary
	arma::ivec dE3;              // RF operator

	// Particle weight:
	arma::vec a;                // Computational particle weigth

	// Initial condition parameters:
	p_IC_TYP p_IC;

	// Boundary conditions:
	p_BC_TYP p_BC;

	// Constructor:
	ionSpecies_TYP(){};

	// Destructor:
	~ionSpecies_TYP(){};
};

//  Define ELECTROMAGNETIC FIELDS DERIVED TYPES:
// =============================================================================
class fields_TYP : public vfield_vec_TYP
{

public:
	vfield_vec_TYP E;   // **** Consider as arma::vec *******
	vfield_vec_TYP B;   // **** Consider as arma::vec *******
  vfield_vec_TYP dB;  // **** Consider as arma::vec *******
  vfield_vec_TYP ddB; // **** Consider as arma::vec *******

	fields_TYP(){};

  /*
	fields_TYP(unsigned int N) : E(N), B(N), dB(N), ddB(N) {};
  */

	~fields_TYP(){};

  /*
	void zeros(unsigned int N);
	void fill(double A);
  */

};

//  Define structure to store simulation geometry information:
// =============================================================================
struct geometry_TYP
{
	double r1;
	double r2;
	double A_0;

	geometry_TYP()
	{
		r1  = 0;
		r2  = 0;
		A_0 = 0;
	}
};

//  Define structure to hold mesh parameters:
// =============================================================================
struct mesh_params_TYP
{
	vfield_vec_TYP nodes;

	int NX_PER_MPI; // Number of mesh nodes along x-axis in subdomain (no ghost nodes considered)
	int NY_PER_MPI; // Number of mesh nodes along y-axis in subdomain (no ghost nodes considered)
	int NZ_PER_MPI; // Number of mesh nodes along z-axis in subdomain (no ghost nodes considered)

	int NX_IN_SIM; // Number of mesh nodes along x-axis in entire simulation domain (no ghost nodes considered)
	int NY_IN_SIM; // Number of mesh nodes along x-axis in entire simulation domain (no ghost nodes considered)
	int NZ_IN_SIM; // Number of mesh nodes along x-axis in entire simulation domain (no ghost nodes considered)

	int NUM_CELLS_IN_SIM; // Number of mesh nodes in the entire simulation domain (no ghost nodes considered)
	int NUM_CELLS_PER_MPI; // Number of mesh nodes in each MPI process (no ghost nodes considered)

	double DX;
	double DY;
	double DZ;

	double LX;		// Size of simulation domain along x-axis
	double LY;		// Size of simulation domain along y-axis
	double LZ;		// Size of simulation domain along z-axis

	// Cartesian unit vectors:
	arma::vec e_x;
	arma::vec e_y;
	arma::vec e_z;

	int SPLIT_DIRECTION;

	mesh_params_TYP()
	{
		e_x = {1.0, 0.0, 0.0};
		e_y = {0.0, 1.0, 0.0};
		e_z = {0.0, 0.0, 1.0};
	}

};

//  Define structure to store fluid species initial condition parameters:
// =============================================================================
struct f_IC_TYP
{
	double ne; 							// Intial reference density
	double Te;
	string Te_fileName;
	int Te_NX;

	// To be used later:
	// ---------------------------------------------------------------------------
	/*
	double Tper;						// Reference perpendicular temperature:
	string Tper_fileName; 			// File containing normalized spatial profile of Tper
	int Tper_NX;						// Number of elements in Tper external file

	double Tpar;						// Reference parallel temperature:
	string Tpar_fileName; 			// File containing normalized spatial profile of Tpar
	int Tpar_NX;						// Number of elements in Tpar external file
	*/
	// ---------------------------------------------------------------------------

	f_IC_TYP()
	{
		ne 			= 0;
		Te      = 0;
		Te_NX   = 0;
	}
};

//  Define structure to store EM field initial condition parameters:
// =============================================================================
struct em_IC_TYP
{
	double BX;
	double BY;
	double BZ;
	string BX_fileName;
	int BX_NX;
	arma::vec Bx_profile;

	double EX;
	double EY;
	double EZ;
	string EX_fileName;
	int EX_NX;
	arma::vec Ex_profile;

	em_IC_TYP()
	{
			BX     = 0;
			BY     = 0;
			BZ     = 0;
			BX_NX  = 0;

			EX     = 0;
			EY     = 0;
			EZ     = 0;
			EX_NX  = 0;
	}
};

//  Define structure to store characteristic values for the normalization:
// =============================================================================
struct CV_TYP
{
	double ne;
	double Te;
	double B;
	double Tpar;
	double Tper;

	CV_TYP()
	{
		ne   = 0;
		Te   = 0;
		B    = 0;
		Tpar = 0;
		Tper = 0;
	}
};

//  Define structure to store switches that control physics modules:
// =============================================================================
struct SW_TYP
{
	int EfieldSolve;
	int HallTermSolve;
	int BfieldSolve;
	int Collisions;
	int RFheating;
	int linearSolve;
	int advancePos;

	SW_TYP()
	{
		EfieldSolve   = 0;
		HallTermSolve = 0;
		BfieldSolve   = 0;
		Collisions    = 0;
		RFheating     = 0;
		linearSolve   = 0;
		advancePos    = 0;
	}

};

//  Define structure to hold MPI parameters:
// =============================================================================
struct mpi_params_TYP
{
	int NUMBER_MPI_DOMAINS;
	int MPI_DOMAIN_NUMBER;
	int FIELDS_ROOT_WORLD_RANK;
	int PARTICLES_ROOT_WORLD_RANK;

	int MPIS_FIELDS;
	int MPIS_PARTICLES;

	int MPI_DOMAINS_ALONG_X_AXIS;
	int MPI_DOMAINS_ALONG_Y_AXIS;
	int MPI_DOMAINS_ALONG_Z_AXIS;

	MPI_Comm MPI_TOPO; // Cartesian topology

	// Particle pusher and field solver communicator params
	int COMM_COLOR;
	int COMM_SIZE;
	int COMM_RANK;
	MPI_Comm COMM;

	bool IS_FIELDS_ROOT;
	bool IS_PARTICLES_ROOT;

	int MPI_CART_COORDS_1D[1];
	int MPI_CART_COORDS_2D[2];
	vector<int *> MPI_CART_COORDS;

	unsigned int iIndex;
	unsigned int fIndex;

	unsigned int irow;
	unsigned int frow;
	unsigned int icol;
	unsigned int fcol;

	int MPI_DOMAIN_NUMBER_CART;
	int LEFT_MPI_DOMAIN_NUMBER_CART;
	int RIGHT_MPI_DOMAIN_NUMBER_CART;
	int UP_MPI_DOMAIN_NUMBER_CART;
	int DOWN_MPI_DOMAIN_NUMBER_CART;
};

// Define structure to hold characteristic scales:
// =============================================================================
struct CS_TYP
{
	double time;
	double velocity;
	double momentum;
	double energy;
	double length;
	double volume;
	double mass;
	double charge;
	double density;
	double eField;
	double bField;
	double pressure;
	double temperature;
	double magneticMoment;
	double resistivity;
	double vacuumPermeability;
	double vacuumPermittivity;

	// Constructor:
	CS_TYP()
	{
		time = 0.0;
		velocity = 0.0;
		momentum = 0.0;
		energy = 0.0;
		length = 0.0;
		volume = 0.0;
		mass = 0.0;
		charge = 0.0;
		density = 0.0;
		eField = 0.0;
		bField = 0.0;
		pressure = 0.0;
		magneticMoment = 0.0;
		vacuumPermeability = 0.0;
		vacuumPermittivity = 0.0;
	}
};

//  Define structure to store simulation parameters:
// =============================================================================
struct params_TYP
{
	// List of variables in the outputs
	vector<string> outputs_variables;

	//int dimensionality;

	//Control parameters for the simulation:
	// Path to save the outputs
	string PATH;

	int argc;
	char **argv;

	// Flag for using a quiet start
	bool quietStart;

	int numberOfRKIterations;
	double smoothingParameter;
	double simulationTime; // In units of the shorter ion gyro-period in the simulation
	int timeIterations;
	double DT;//Time step
	double DTc;//Cyclotron period fraction.

	// Consider deleting of not neded:
	// ===============================
	//int loadFields;
	int loadGrid;
	int usingHDF5;
	int outputCadenceIterations;
	arma::file_type outputFormat;//Outputs format (raw_ascii,raw_binary).
	//types_info typesInfo;
	// ===============================

	double outputCadence;//Save variables each "outputCadence" times the background ion cycloperiod.

	// Geometric information:
	geometry_TYP geometry;

	// Mesh parameters:
	mesh_params_TYP mesh;

	// Ions properties:
	int numberOfParticleSpecies; // This species are evolved self-consistently with the fields
	int numberOfTracerSpecies; // This species are not self-consistently evolved with the fields

	// Electron initial conditions:
	f_IC_TYP f_IC;

	// EM initial conditions:
	em_IC_TYP em_IC;

	// Simulation Characterstic values:
	CV_TYP CV;

	// Simulation switches:
	SW_TYP SW;

	int filtersPerIterationFields;
	int filtersPerIterationIons;

	double ionLarmorRadius;
	double ionSkinDepth;
	double ionGyroPeriod;

	double DrL;
	double dp;

	int checkStability;
	int rateOfChecking;//Each 'rateOfChecking' iterations we use the CLF criteria for the particles to stabilize the simulation

	// MPI parameters
	mpi_params_TYP mpi;

	// Error codes
	map<int,string> errorCodes;

	// Constructor
	params_TYP(){};
};

// Define structure to hold fundamental scales:
// =============================================================================
struct FS_TYP
{
	double electronSkinDepth;
	double electronGyroPeriod;
	double electronGyroRadius;
	double * ionSkinDepth;
	double * ionGyroPeriod;
	double * ionGyroRadius;

	// Constructor:
	FS_TYP(params_TYP * params)
	{
		electronSkinDepth = 0.0;
		electronGyroPeriod = 0.0;
		electronGyroRadius = 0.0;
		ionSkinDepth = new double[params->numberOfParticleSpecies];
		ionGyroPeriod = new double[params->numberOfParticleSpecies];
		ionGyroRadius = new double[params->numberOfParticleSpecies];
	}
};

#endif