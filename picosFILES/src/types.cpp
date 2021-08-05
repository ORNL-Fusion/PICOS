#include "types.h"

double F_EPSILON_DS = F_EPSILON;    // Dimensionless vacuum permittivity
double F_E_DS = F_E;                // Dimensionless electric charge
double F_ME_DS = F_ME;              // Dimensionless electron mass
double F_MU_DS = F_MU;              // Dimensionless vacuum permeability
double F_C_DS = F_C;                // Dimensionless speed of light

// fields_TYP functions:
// =============================================================================
void fields_TYP::zeros(unsigned int N)
{
	EX_m.zeros(N);
	BX_m.zeros(N);
	dBX_m.zeros(N);
	ddBX_m.zeros(N);
}

void fields_TYP::fill(double A)
{
	EX_m.fill(A);
	BX_m.fill(A);
	dBX_m.fill(A);
	ddBX_m.fill(A);
}

/*
// vfield_vec_TYP functions:
// =============================================================================
// Constructors:
vfield_vec_TYP::vfield_vec_TYP(unsigned int N)
{
	X = arma::zeros(N);
	Y = arma::zeros(N);
	Z = arma::zeros(N);
}

void vfield_vec_TYP::ones(unsigned int N)
{
	X.ones(N);
	Y.ones(N);
	Z.ones(N);
}

void vfield_vec_TYP::zeros()
{
	X.zeros();
	Y.zeros();
	Z.zeros();
}

void vfield_vec_TYP::zeros(unsigned int N){
	X.zeros(N);
	Y.zeros(N);
	Z.zeros(N);
}

void vfield_vec_TYP::fill(double value)
{
	X.fill(value);
	Y.fill(value);
	Z.fill(value);
}

*/
