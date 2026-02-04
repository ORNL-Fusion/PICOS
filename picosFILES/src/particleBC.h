#ifndef H_PARTICLEBOUNDARYCONDITIONS
#define H_PARTICLEBOUNDARYCONDITIONS

#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <cassert>
#include <numbers>

#define ARMA_ALLOW_FAKE_GCC
#include "armadillo"
#include "types.h"
#include "mpi_main.h"

#include "parallel_random.hpp"

using namespace std;
using namespace arma;

class particleBC_TYP
{

private:
    typedef uniform_real_distribution<double> uniform;
    typedef picos::random::instance<double, uniform, 0.0, 2*numbers::pi_v<double>> uniform_2Pi;
    typedef picos::random::instance<double, uniform, 0.0, 1.0> uniform_one;
    std::random_device device;
    std::vector<uniform_2Pi> randoms_2pi;
    std::vector<uniform_one> randoms_one;

    void particleReinjection(int ii, const params_TYP &params, const CS_TYP &CS, fields_TYP &fields, ionSpecies_TYP &ION, uniform_2Pi &rand_2pi, uniform_one &rand_one) const;

    template <size_t S>
    static void MPI_AllreduceDouble(const params_TYP &params, double *v)
    {
        MPI_Allreduce(MPI_IN_PLACE, v, S, MPI_DOUBLE, MPI_SUM, params.mpi.COMM);
    }

    template <typename vec_TYP>
    static void MPI_OMP_AllreduceVec(const params_TYP &params, vec_TYP &V1, vec_TYP &V2, std::array<double,2> &S)
    {
        // Clear S:
        // =======
        S = {0,0};

        // Reduce S over all threads in a single MPI process:
        // ==================================================
        // NOTE: Assumes v1 and v2 are the same size.
        assert(V1.n_elem == V2.n_elem && "V1 and V2 have different sizes.");
        #pragma omp declare reduction(sum: std::array<double, 2> : omp_out[0] += omp_in[0], omp_out[1] += omp_in[1])
        #pragma omp parallel for default(none) shared(V, std::cout) reduction(sum:S)
        for(int ii=0, iie=V1.n_elem; ii<iie; ii++)
        {
            S[0] += V1(ii);
            S[1] += V2(ii);
        }

        // AllReduce S over all PARTICLE MPIs:
        // ==================================
        MPI_AllreduceDouble<2> (params,S.data());
    }

    void checkBoundaryAndFlag(const params_TYP &params, const CS_TYP &CS, fields_TYP &fields, vector<ionSpecies_TYP> &IONS) const;

    void getFluxesAcrossBoundaries(const params_TYP &params, const CS_TYP &CS, fields_TYP &fields, vector<ionSpecies_TYP> &IONS);

    void calculateParticleWeight(const params_TYP &params, const CS_TYP &CS, fields_TYP &fields, vector<ionSpecies_TYP> &IONS) const;

    void getParticleInjectionRates(const params_TYP &params, const CS_TYP &CS, fields_TYP &fields, vector<ionSpecies_TYP> &ION);

public:

    struct dot_buffer {
        double N1;
        double E1;
        double N2;
        double E2;
        double N5;
        double E5;
    };
    dot_buffer dot_;

    particleBC_TYP();

    void applyParticleReinjection(const params_TYP &params, const CS_TYP &CS, fields_TYP &fields, vector<ionSpecies_TYP> &IONS);
};

#endif
