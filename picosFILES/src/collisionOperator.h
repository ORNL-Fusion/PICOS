// FP_TYP::ApplyCollisions_AllSpecies()
#ifndef H_COLLISIONOPERATOR
#define H_COLLISIONOPERATOR

#include <iostream>
#include <cmath>
#include <vector>

#define ARMA_ALLOW_FAKE_GCC
#include "armadillo"
#include "types.h"
#include "mpi_main.h"

#include "parallel_random.hpp"

using namespace std;
using namespace arma;

class coll_operator_TYP
{
private:
    // A uniform distribution.
    typedef uniform_int_distribution<short> uniform;
    // A Uniform random instance type.
    typedef picos::random::instance<short, uniform, 0, 1> uniform_random;

    // Ion moment interpolation functions:
    void interpolateIonMoments(const params_TYP * params, vector<ionSpecies_TYP> * IONS, int a, int b);
    void interpolateScalarField(const params_TYP * params, ionSpecies_TYP * IONS, arma::vec * F_m, arma::vec * F_p);
    void fill4Ghosts(arma::vec * v);

    // Electron temperature interpolation:
    void interpolateElectronTemperature(const params_TYP * params, vector<ionSpecies_TYP> * IONS, int a, electrons_TYP * electrons);

    // Scattering operators:
    void u_CollisionOperator(double &w, const double xab, const double wTb,
                             const double nb, const double Tb,
                             const double Mb, const double Zb,
                             const double Za, const double Ma,
                             const double DT, uniform_random &rand);
    void xi_CollisionOperator(double &xi, const double xab, const double wTb,
                              const double nb, const double Tb,
                              const double Mb, const double Zb,
                              const double Za, const double Ma,
                              const double DT, uniform_random &rand);

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

    std::random_device device;
    std::vector<uniform_random> randoms;

public:
    coll_operator_TYP() :
    randoms(picos::random::instances<short, uniform, 0, 1> (device())) {}

    void ApplyCollisions_AllSpecies(const params_TYP * params, const CS_TYP * CS, vector<ionSpecies_TYP> * IONS, electrons_TYP * electrons);
};

#endif
