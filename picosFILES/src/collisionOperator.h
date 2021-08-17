// FP_TYP::ApplyCollisions_AllSpecies()
#ifndef H_COLLISIONOPERATOR
#define H_COLLISIONOPERATOR

#include <iostream>
#include <cmath>
#include <vector>

#include "armadillo"
#include "types.h"
#include "mpi_main.h"
#include "omp.h"

using namespace std;
using namespace arma;

class coll_operator_TYP
{
private:
    // Ion moment interpolation functions:
    void interpolateIonMoments(const params_TYP * params, vector<ionSpecies_TYP> * IONS, int a, int b);
    void interpolateScalarField(const params_TYP * params, ionSpecies_TYP * IONS, arma::vec * F_m, arma::vec * F_p);
    void fill4Ghosts(arma::vec * v);

    // Scattering operators:
    void u_CollisionOperator(double * w, double xab, double wTb, double nb, double Tb, double Mb, double Zb, double Za, double Ma, double DT);
    void xi_CollisionOperator(double * xi, double xab, double wTb, double nb, double Tb, double Mb, double Zb, double Za, double Ma, double DT);

    /*
    void phi_CollisionOperator(double * phi, double xi, double xab, double wTb, double nb, double Tb, double Mb, double Zb, double Za, double Ma, double DT);
    */

    // Coordinate transformation:
    void cartesian2Spherical(double * wx, double * wy, double * wz, double * w, double * xi, double * phi);
    void Spherical2Cartesian(double * w, double * xi, double * phi, double * wx, double * wy, double * wz);

    // Coulomb colliional rates:
    double nu_E(double xab, double nb, double Tb, double Mb, double Zb, double Za, double Ma, int energyOperatorModel);
    double nu_D(double xab, double nb, double Tb, double Mb, double Zb, double Za, double Ma);
    double nu_ab0(double nb, double Tb, double Mb, double Zb, double Za, double Ma);
    double logA(double nb, double Tb);
    double Gb(double xab);
    double erfp(double xab);
    double erfpp(double xab);
    double E_nuE_d_nu_E_dE(double xab);

public:
    coll_operator_TYP();
    void ApplyCollisions_AllSpecies(const params_TYP * params, const CS_TYP * CS, vector<ionSpecies_TYP> * IONS);
};

#endif
