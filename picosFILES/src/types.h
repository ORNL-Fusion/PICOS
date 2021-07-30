#ifndef H_TYPES
#define H_TYPES

#include <vector>
#include <armadillo>

// Declare vector field types:
// =============================================================================
class vfield_vec_TYP
{

public:
	arma::vec X;
	arma::vec Y;
	arma::vec Z;

	vfield_vec_TYP(){};

  /*
	vfield_vec_TYP(unsigned int N);
  */

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
	std::string Tper_fileName; 				// File containing normalized spatial profile of Tper
	std::string Tpar_fileName; 				// File containing normalized spatial profile of Tpar
	std::string densityFraction_fileName;	// File containing normalized spatial profile

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
	std::string G_fileName;

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

  /*
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

	double densityFraction;		//

	double go;					// Initial relativistic gamma
	double LarmorRadius;		// Larmor radius.
	double VTper;				// Thermal velocity.
	double VTpar;				// Thermal velocity.
	double Wc;					// Average cyclotron frequency.
	double Wp;					// Plasma frequency.
	double avg_mu; 				// Average magnetic moment
  */

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

#endif
