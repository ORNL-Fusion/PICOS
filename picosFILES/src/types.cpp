#include "types.h"

double F_EPSILON_DS = F_EPSILON;    // Dimensionless vacuum permittivity
double F_E_DS = F_E;                // Dimensionless electric charge
double F_ME_DS = F_ME;              // Dimensionless electron mass
double F_MU_DS = F_MU;              // Dimensionless vacuum permeability
double F_C_DS = F_C;                // Dimensionless speed of light

// Constructors:
vfield_vec_TYP::vfield_vec_TYP(unsigned int N)
{
	X = arma::zeros(N);
	Y = arma::zeros(N);
	Z = arma::zeros(N);
}
