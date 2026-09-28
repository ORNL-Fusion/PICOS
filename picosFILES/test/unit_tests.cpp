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
    params.pairSource.weightMode = PAIR_SOURCE_WEIGHT_EXPLICIT_RATE;
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
    allocate_test_species(species[1], -1.0, 40.0);
    species[0].X_p(0) = -0.1;
    species[0].X_p(1) = 0.5;
    species[1].X_p(0) = 1.1;
    species[1].X_p(1) = 0.5;

    particleBC_TYP bc;
    bc.applyParticleReinjection(params, cs, fields, species);

    const double expectedIonWeight = params.pairSource.rate*params.DT/species[0].NCP;
    const double expectedElectronWeight = params.pairSource.rate*params.DT/species[1].NCP;
    assert(std::abs(species[0].X_p(0) - params.pairSource.mean_x) < 1.0e-15);
    assert(std::abs(species[1].X_p(0) - params.pairSource.mean_x) < 1.0e-15);
    assert(std::abs(species[0].X_p(0) - species[1].X_p(0)) < 1.0e-15);
    assert(std::abs(species[0].a_p(0) - expectedIonWeight) < 1.0e-15);
    assert(std::abs(species[1].a_p(0) - expectedElectronWeight) < 1.0e-15);
    assert(std::abs(species[0].Z*species[0].NCP*species[0].a_p(0) +
                    species[1].Z*species[1].NCP*species[1].a_p(0)) < 1.0e-15);
    assert(species[0].f1(0) == 0 && species[0].f2(0) == 0);
    assert(species[1].f1(0) == 0 && species[1].f2(0) == 0);

    // With a fixed physical source rate, halving DT must halve the number of
    // physical pairs injected while preserving charge balance for unequal NCP.
    params.DT *= 0.5;
    std::vector<ionSpecies_TYP> half_step_species(2);
    allocate_test_species(half_step_species[0], 1.0, 10.0);
    allocate_test_species(half_step_species[1], -1.0, 40.0);
    half_step_species[0].X_p(0) = -0.1;
    half_step_species[0].X_p(1) = 0.5;
    half_step_species[1].X_p(0) = 1.1;
    half_step_species[1].X_p(1) = 0.5;
    particleBC_TYP half_step_bc;
    half_step_bc.applyParticleReinjection(params, cs, fields, half_step_species);
    const double full_step_pairs = species[0].NCP*species[0].a_p(0);
    const double half_step_ion_pairs = half_step_species[0].NCP*half_step_species[0].a_p(0);
    const double half_step_electron_pairs = half_step_species[1].NCP*half_step_species[1].a_p(0);
    assert(std::abs(half_step_ion_pairs - 0.5*full_step_pairs) < 1.0e-15);
    assert(std::abs(half_step_electron_pairs - half_step_ion_pairs) < 1.0e-15);
}

void reformulated_logical_sheath_unit_test() {
    params_TYP params;
    params.mpi.COMM = MPI_COMM_WORLD;
    params.mpi.COMM_COLOR = PARTICLES_MPI_COLOR;
    params.mpi.IS_PARTICLES_ROOT = true;
    params.SW.fieldSolveModel = FIELD_SOLVE_REFORMULATED_POISSON;
    params.SW.pairSource = 0;
    params.em_IC.poissonBCModel = POISSON_BC_SHEATH;
    params.em_IC.sheathCoefficient = 3.0;
    params.f_IC.Te = F_E_DS;
    params.advanceParticleMethod = PARTICLE_PUSH_GC_VPER;
    params.DT = 0.01;
    params.geometry.LX_min = 0.0;
    params.geometry.LX_max = 1.0;
    params.geometry.LX = 1.0;
    params.mesh.DX = 0.1;
    params.mesh.NX_IN_SIM = 10;

    fields_TYP fields;
    fields.Phi_m.zeros(params.mesh.NX_IN_SIM + 2);
    CS_TYP cs;
    std::vector<ionSpecies_TYP> species(1);
    allocate_test_species(species[0], -1.0, 1.0);
    species[0].X_p(0) = -0.01;
    species[0].X_p(1) = 0.5;
    species[0].V_p(0,0) = -0.1;

    particleBC_TYP bc;
    bc.applyParticleReinjection(params, cs, fields, species);

    assert(std::abs(species[0].X_p(0) - 0.25*params.mesh.DX) < 1.0e-15);
    assert(std::abs(species[0].V_p(0,0) - 0.1) < 1.0e-15);
    assert(species[0].f1(0) == 0 && species[0].f2(0) == 0);
    assert(std::abs(species[0].a_p(0) - 1.0) < 1.0e-15);
}

void current_balanced_logical_sheath_unit_test() {
    params_TYP params;
    params.mpi.COMM = MPI_COMM_WORLD;
    MPI_Comm_size(MPI_COMM_WORLD, &params.mpi.COMM_SIZE);
    params.mpi.COMM_COLOR = PARTICLES_MPI_COLOR;
    params.mpi.IS_PARTICLES_ROOT = true;
    params.SW.fieldSolveModel = FIELD_SOLVE_REFORMULATED_POISSON;
    params.SW.pairSource = 0;
    params.em_IC.poissonBCModel = POISSON_BC_SHEATH;
    params.em_IC.sheathCurrentBalance = 1;
    params.advanceParticleMethod = PARTICLE_PUSH_GC_VPER;
    params.DT = 0.01;
    params.geometry.LX_min = 0.0;
    params.geometry.LX_max = 1.0;
    params.geometry.LX = 1.0;
    params.mesh.DX = 0.1;
    params.mesh.NX_IN_SIM = 10;

    fields_TYP fields;
    fields.Phi_m.zeros(params.mesh.NX_IN_SIM + 2);
    CS_TYP cs;
    std::vector<ionSpecies_TYP> species(2);
    allocate_test_species(species[0], 1.0, 2.0);
    allocate_test_species(species[1], -1.0, 2.0);

    // One ion carries two units of outgoing charge.  Of the two equal-weight
    // electrons, only the higher-energy one should pass the logical sheath.
    species[0].X_p(0) = -0.01;
    species[0].X_p(1) = 0.5;
    species[1].X_p(0) = -0.01;
    species[1].X_p(1) = -0.01;
    species[1].V_p(0,0) = -0.1;
    species[1].V_p(1,0) = -1.0;

    particleBC_TYP bc;
    bc.applyParticleReinjection(params, cs, fields, species);

    assert(std::abs(species[1].X_p(0) - 0.25*params.mesh.DX) < 1.0e-15);
    assert(species[1].V_p(0,0) > 0.0);
    assert(std::abs(species[1].X_p(1) - 0.25*params.mesh.DX) > 1.0e-12);

    // Verify that ion charge is carried across particle-sparse steps.  First
    // lose an ion with no electron candidate, then present one equal-charge
    // electron on the following step; it must be transmitted, not reflected.
    std::vector<ionSpecies_TYP> sparseSpecies(2);
    allocate_test_species(sparseSpecies[0], 1.0, 2.0);
    allocate_test_species(sparseSpecies[1], -1.0, 2.0);
    sparseSpecies[0].X_p(0) = -0.01;
    sparseSpecies[0].X_p(1) = 0.5;
    sparseSpecies[1].X_p.fill(0.5);
    particleBC_TYP sparseBc;
    sparseBc.applyParticleReinjection(params, cs, fields, sparseSpecies);

    sparseSpecies[0].X_p.fill(0.5);
    sparseSpecies[1].X_p(0) = -0.01;
    sparseSpecies[1].X_p(1) = 0.5;
    sparseSpecies[1].V_p(0,0) = -0.1;
    sparseBc.applyParticleReinjection(params, cs, fields, sparseSpecies);
    assert(std::abs(sparseSpecies[1].X_p(0) - 0.25*params.mesh.DX) > 1.0e-12);
}

void sonic_bohm_outflow_unit_test() {
    params_TYP params;
    params.mpi.COMM = MPI_COMM_WORLD;
    params.mpi.COMM_COLOR = PARTICLES_MPI_COLOR;
    params.numberOfParticleSpecies = 1;
    params.SW.Bohm = 1;
    params.bohm.type = 2;
    params.bohm.edgeCells = 1;
    params.bohm.tOn = 0.0;
    params.bohm.gammaI = 3.0;
    params.currentTime = 0.0;
    params.geometry.LX_min = 0.0;
    params.geometry.LX_max = 1.0;
    params.mesh.DX = 0.1;
    params.mesh.NX_IN_SIM = 10;

    CS_TYP cs;
    cs.time = 1.0;
    electrons_TYP electrons;
    electrons.Te_m.ones(12);
    electrons.Te_m *= 2.0;

    std::vector<ionSpecies_TYP> species(1);
    allocate_test_species(species[0], 1.0, 1.0);
    species[0].M = 2.0;
    species[0].X_p(0) = 0.05;
    species[0].X_p(1) = 0.95;
    species[0].V_p(0,0) = -0.25; // subsonic left outflow
    species[0].V_p(1,0) = 2.0;   // supersonic right outflow

    particleBC_TYP bc;
    bc.enforceSonicBohmOutflow(params, cs, electrons, species);

    // With Te/Mi = 1 and zero parallel thermal spread, cs = 1.
    assert(std::abs(species[0].V_p(0,0) + 1.0) < 1.0e-14);
    assert(std::abs(species[0].V_p(1,0) - 2.0) < 1.0e-14);

    // Reverse which side is subsonic: the supersonic left side must remain
    // unchanged and the subsonic right side must be shifted to +cs.
    species[0].V_p(0,0) = -2.0;
    species[0].V_p(1,0) = 0.25;
    bc.enforceSonicBohmOutflow(params, cs, electrons, species);
    assert(std::abs(species[0].V_p(0,0) + 2.0) < 1.0e-14);
    assert(std::abs(species[0].V_p(1,0) - 1.0) < 1.0e-14);
}

void kinetic_electron_thermostat_unit_test() {
    params_TYP params;
    params.mpi.COMM = MPI_COMM_WORLD;
    params.mpi.COMM_COLOR = PARTICLES_MPI_COLOR;
    params.SW.kineticElectronThermostat = 1;
    params.kineticElectronThermostatRelaxation = 1.0;
    params.advanceParticleMethod = PARTICLE_PUSH_GC_VPER;
    params.mesh.NX_IN_SIM = 1;

    std::vector<ionSpecies_TYP> species(1);
    allocate_test_species(species[0], -1.0, 1.0);
    species[0].M = 1.0;
    species[0].mn.zeros(2);
    species[0].Te_p.ones(2);
    species[0].Te_p *= 4.0;
    species[0].BX_p.ones(2);
    species[0].mu_p.zeros(2);
    species[0].V_p(0,0) = 0.0;
    species[0].V_p(1,0) = 2.0;
    species[0].V_p.col(1).fill(2.0);

    const double thermostatPower = PIC_TYP::applyKineticElectronThermostat(params, species);

    const double meanVpar = arma::mean(species[0].V_p.col(0));
    const double Tpar = arma::mean(arma::square(species[0].V_p.col(0) - meanVpar));
    const double Tper = 0.5*arma::mean(arma::square(species[0].V_p.col(1)));
    assert(std::abs(meanVpar - 1.0) < 1.0e-14);
    assert(std::abs(Tpar - 4.0) < 1.0e-14);
    assert(std::abs(Tper - 4.0) < 1.0e-14);
    assert(thermostatPower > 0.0);
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
        reformulated_logical_sheath_unit_test();
        current_balanced_logical_sheath_unit_test();
        sonic_bohm_outflow_unit_test();
        kinetic_electron_thermostat_unit_test();
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
