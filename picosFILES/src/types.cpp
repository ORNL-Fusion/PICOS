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
