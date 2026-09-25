//------------------------------------------------------------------------------
///  @file unit_test.cpp
///  @brief Runs the picos unit tests.
//------------------------------------------------------------------------------

//  Turn on asserts even in release builds.
#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../src/collisionOperator.h"
#include "../src/particleBC.h"
#include "../src/PIC.h"

namespace {
void guiding_center_mirror_force_unit_test() {
    const double qa = 2.0;
    const double ma = 4.0;
    const std::array<double, 3> fields = {3.0, 5.0, 7.0};
    const std::array<double, 3> state = {0.0, 11.0, 13.0};
    std::array<double, 3> rhs{};

    PIC_TYP::guidingCenterVperRhs(qa, ma, fields, state, rhs);

    const double expectedParallel =
        -0.5*state[2]*state[2]*fields[2]/fields[1] + qa*fields[0]/ma;
    assert(std::abs(rhs[0] - state[1]) < 1.0e-15);
    assert(std::abs(rhs[1] - expectedParallel) < 1.0e-13);
    assert(std::abs(rhs[2] - 0.5*state[2]*state[1]*fields[2]/fields[1]) < 1.0e-13);

    // The parallel mirror force depends on v_perpendicular squared and is
    // therefore unchanged when v_parallel reverses direction.
    auto reversed = state;
    reversed[1] *= -1.0;
    PIC_TYP::guidingCenterVperRhs(qa, ma, fields, reversed, rhs);
    assert(std::abs(rhs[1] - expectedParallel) < 1.0e-13);
}

void allocate_test_species(ionSpecies_TYP& species, double z, double ncp) {
    species.Z = z;
    species.Q = z;
    species.M = z < 0.0 ? F_ME : 2.0*F_MP;
    species.NSP = 2;
    species.NCP = ncp;
    species.X_p.zeros(2);
    species.V_p.zeros(2, 2);
    species.a_p.ones(2);
    species.f1.zeros(2);
    species.f2.zeros(2);
    species.f3.zeros(2);
    species.f5.zeros(2);
    species.dE1.zeros(2);
    species.dE2.zeros(2);
    species.dE3.zeros(2);
    species.dE5.zeros(2);
    species.p_BC.BC_type = 1;
    species.p_BC.T = 1.0;
    species.p_BC.E = 0.0;
    species.p_BC.mean_x = 0.5;
    species.p_BC.sigma_x = 0.1;
}

void pair_source_unit_test() {
    params_TYP params;
    params.mpi.COMM = MPI_COMM_WORLD;
    params.mpi.COMM_COLOR = PARTICLES_MPI_COLOR;
    params.mpi.IS_PARTICLES_ROOT = true;
    params.SW.fieldSolveModel = FIELD_SOLVE_OHM;
    params.SW.pairSource = 1;
    params.advanceParticleMethod = PARTICLE_PUSH_GC_VPER;
    params.DT = 0.01;
    params.geometry.LX_min = 0.0;
    params.geometry.LX_max = 1.0;
    params.geometry.LX = 1.0;
    params.mesh.DX = 0.1;
    params.pairSource.ionSpecies = 0;
    params.pairSource.electronSpecies = 1;
    params.pairSource.rate = 100.0;
    params.pairSource.mean_x = 0.5;
    params.pairSource.sigma_x = 0.0;
    params.pairSource.ionT = 1.0;
    params.pairSource.electronT = 1.0;
    params.pairSource.ionE = 0.0;
    params.pairSource.electronE = 0.0;
    params.pairSource.ionEta = 0.0;
    params.pairSource.electronEta = 0.0;
    params.pairSource.maxParticleWeight = 1000.0;

    fields_TYP fields;
    CS_TYP cs;
    std::vector<ionSpecies_TYP> species(2);
    allocate_test_species(species[0], 1.0, 10.0);
    allocate_test_species(species[1], -1.0, 10.0);
    species[0].X_p(0) = -0.1;
    species[0].X_p(1) = 0.5;
    species[1].X_p(0) = 1.1;
    species[1].X_p(1) = 0.5;

    particleBC_TYP bc;
    bc.applyParticleReinjection(params, cs, fields, species);

    const double expectedWeight = params.pairSource.rate*params.DT/species[0].NCP;
    assert(std::abs(species[0].X_p(0) - params.pairSource.mean_x) < 1.0e-15);
    assert(std::abs(species[1].X_p(0) - params.pairSource.mean_x) < 1.0e-15);
    assert(std::abs(species[0].X_p(0) - species[1].X_p(0)) < 1.0e-15);
    assert(std::abs(species[0].a_p(0) - expectedWeight) < 1.0e-15);
    assert(std::abs(species[1].a_p(0) - expectedWeight) < 1.0e-15);
    assert(std::abs(species[0].Z*species[0].NCP*species[0].a_p(0) +
                    species[1].Z*species[1].NCP*species[1].a_p(0)) < 1.0e-15);
    assert(species[0].f1(0) == 0 && species[0].f2(0) == 0);
    assert(species[1].f1(0) == 0 && species[1].f2(0) == 0);
}
}

//------------------------------------------------------------------------------
///  @brief Run tests.
///
///  @tparam T Base type of the calculation.
//------------------------------------------------------------------------------
template<std::floating_point T> void run_tests() {
    picos::random::test<T> ();

    if constexpr (std::is_same<T, double> ()) {
        coll_operator_TYP opt;
        opt.unit_test();
        pair_source_unit_test();
        guiding_center_mirror_force_unit_test();
    }
}

//------------------------------------------------------------------------------
///  @brief Main program of the test.
///
///  @param[in] argc Number of commandline arguments.
///  @param[in] argv Array of commandline arguments.
//------------------------------------------------------------------------------
int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);

    run_tests<float> ();
    run_tests<double> ();

    MPI_Finalize();
}
