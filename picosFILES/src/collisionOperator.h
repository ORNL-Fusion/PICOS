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
    void interpolateIonMoments(const params_TYP &params, ionSpecies_TYP &iona, const ionSpecies_TYP &ionb) const;
    void interpolateScalarField(const params_TYP &params, const ionSpecies_TYP &ion, const arma::vec &F_m, arma::vec &F_p) const;
    void fill4Ghosts(arma::vec &v) const;

    // Electron temperature interpolation:
    void interpolateElectronTemperature(const params_TYP &params, ionSpecies_TYP &ion, const electrons_TYP &electrons) const;

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

    // Collisional rates based on Maxwellian background species:
    // =============================================================================
    template<uint8_t energyOperatorModel=2>
    double nu_E(const double xab, const double nb, const double Tb, const double Mb, const double Zb, const double Za, const double Ma) const
    {
        const double mass_ratio = 2.0*Ma/Mb*Gb(xab);
        const double nu = nu_ab0(nb,Tb,Mb,Zb,Za,Ma)/xab;
        if constexpr (energyOperatorModel == 1)
        {
            // From Hinton 1983 EQ 92 and T.S. Chen 1988 EQ 50
            return nu*(mass_ratio - erfp(xab)/xab);
        }
        else if constexpr (energyOperatorModel == 2)
        {
            //From T.S. Chen 1988 Report EQ 57 commonly used for NBI
            return nu*mass_ratio;
        }
        static_assert(energyOperatorModel != 1 ||
                      energyOperatorModel != 2,
                      "Invalid energy operator model.");

        /* References:
        T.S Chen 1988:
        "A General Form of the Coulomb Scattering Operators for Monte Carlo ...
        Simulations and a Note on the Guiding Center Equations in Different Magnetic Coordinate Conventions"

        Hinton 1983:
        "Handbook of Plasma Physics
        Editors: M.N. Rosenbluth and R.Z. Sagdeev
        Chapter 1.5 - Collisional Transport in Plasma"
        */
    }

    double nu_D(const double xab, const double nb, const double Tb, const double Mb, const double Zb, const double Za, const double Ma) const;
    double nu_ab0(const double nb, const double Tb, const double Mb, const double Zb, const double Za, const double Ma) const;
    double logA(const double nb, const double Tb) const;
    double Gb(const double xab) const;
    double erfp(const double xab) const;
    double erfpp(double xab) const;
    double E_nuE_d_nu_E_dE(const double xab) const;

    std::random_device device;
    std::vector<uniform_random> randoms;

public:
    coll_operator_TYP() :
    randoms(picos::random::instances<short, uniform, 0, 1> (device())) {}

    void ApplyCollisions_AllSpecies(const params_TYP &params, const CS_TYP &CS, vector<ionSpecies_TYP> &IONS, electrons_TYP &electrons);
};

#endif
