#!/usr/bin/env python3
"""Build and plot the MPEX scenario-14 Fig. 17-style E-z distribution.

The figure quantity is the electron distribution in axial position and kinetic
energy, with the scenario-14 magnetic-field profile overlaid on a secondary
axis.  This script uses the checked-in MPEX scenario-14 template profiles and
can optionally run the Fortran reference and PICOS++ decks it generates.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
from pathlib import Path
from typing import Any

import numpy as np

from compare_fortran_picos_short import (
    C_LIGHT,
    E_CHARGE,
    M_E,
    ensure_plotting,
    parse_picos_input,
    picos_components,
    read_fortran_record_float64,
    read_hdf5_dataset,
    read_hdf5_scalar_optional,
)


def ion_gyroperiod(b_t: float, mass_amu: float = 2.0, charge_state: float = 1.0) -> float:
    atomic_mass = 1.66053906660e-27
    electron_charge = 1.602176634e-19
    mass_kg = mass_amu * atomic_mass
    omega_ci = abs(charge_state) * electron_charge * b_t / mass_kg
    return 2.0 * math.pi / omega_ci


def ion_skin_depth(ne_m3: float, mass_amu: float = 2.0, charge_state: float = 1.0) -> float:
    atomic_mass = 1.66053906660e-27
    electron_charge = 1.602176634e-19
    eps0 = 8.8541878128e-12
    c_light = 299792458.0
    mass_kg = mass_amu * atomic_mass
    omega_pi = math.sqrt(ne_m3 * (charge_state * electron_charge) ** 2 / (eps0 * mass_kg))
    return c_light / omega_pi


def scenario14_profile(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    raw = np.loadtxt(args.scenario14_b_file)
    z = np.linspace(args.z_min, args.z_max, raw.size)
    b_t = raw * args.b_scale
    b_norm = b_t / np.max(np.abs(b_t))
    return z, b_t, b_norm


def rf_frequency_for_resonance_b(b_res_t: float, n_harmonic: int) -> float:
    electron_charge = 1.602176634e-19
    electron_mass = 9.1093837015e-31
    return electron_charge * b_res_t * n_harmonic / (2.0 * math.pi * electron_mass)


def resonance_b(args: argparse.Namespace) -> float:
    if args.resonance_b_t is not None:
        return float(args.resonance_b_t)
    electron_charge = 1.602176634e-19
    electron_mass = 9.1093837015e-31
    return 2.0 * math.pi * electron_mass * args.rf_frequency / (electron_charge * args.n_harmonic)


def resonance_crossings(z: np.ndarray, b_t: np.ndarray, b_res: float) -> list[float]:
    crossings: list[float] = []
    for ii in range(z.size - 1):
        y0 = b_t[ii] - b_res
        y1 = b_t[ii + 1] - b_res
        if y0 == 0.0:
            crossings.append(float(z[ii]))
        if y0 * y1 < 0.0:
            crossings.append(float(z[ii] + (0.0 - y0) * (z[ii + 1] - z[ii]) / (y1 - y0)))
    return crossings


def selected_resonance_z(args: argparse.Namespace, z: np.ndarray, b_t: np.ndarray) -> float:
    if args.resonance_z is not None:
        return float(args.resonance_z)
    b_res = resonance_b(args)
    crossings = resonance_crossings(z, b_t, b_res)
    if not crossings:
        return float(z[np.argmin(np.abs(b_t - b_res))])
    return min(crossings, key=lambda value: abs(value - args.reference_resonance_z))


def write_vector(path: Path, values: np.ndarray) -> None:
    np.savetxt(path, values, fmt="%.16e")


def picos_profile_z(args: argparse.Namespace) -> np.ndarray:
    """Coordinate used by PICOS++ for external B-field profile files.

    initializeElectromagneticFields() interprets the external B file as two
    guard-like points wider than the active mesh:
    z_i = LX_min - 0.5*dz + i*dz, dz = LX/(Nfile - 2).
    """
    dx = (args.z_max - args.z_min) / (args.profile_points - 2)
    return args.z_min - 0.5 * dx + dx * np.arange(args.profile_points)


def picos_aux_profile_z(args: argparse.Namespace) -> np.ndarray:
    """Coordinate used by PICOS++ for IC and pair-source profile files."""
    return np.linspace(args.z_min, args.z_max, args.profile_points)


def normalize_shape(values: np.ndarray, floor_fraction: float = 0.0) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    values = np.where(np.isfinite(values), values, 0.0)
    values = np.maximum(values, 0.0)
    vmax = float(values.max()) if values.size else 0.0
    if vmax <= 0.0:
        return np.ones_like(values)
    shape = values / vmax
    if floor_fraction > 0.0:
        shape = np.maximum(shape, floor_fraction)
    return shape


def gaussian_shape(z: np.ndarray, center: float, sigma: float, floor_fraction: float = 0.0) -> np.ndarray:
    sigma = max(float(sigma), 1.0e-6)
    shape = np.exp(-0.5 * np.square((z - center) / sigma))
    if floor_fraction > 0.0:
        shape = np.maximum(shape, floor_fraction)
    return normalize_shape(shape)


def read_csv_column(data: np.ndarray, *names: str) -> np.ndarray | None:
    if data.dtype.names is None:
        return None
    names_by_lower = {name.lower(): name for name in data.dtype.names}
    for name in names:
        actual = names_by_lower.get(name.lower())
        if actual is not None:
            return np.asarray(data[actual], dtype=float)
    return None


def profile_from_particles_nc(path: Path, z_profile: np.ndarray) -> np.ndarray | None:
    try:
        import netCDF4 as nc
    except Exception:
        return None

    with nc.Dataset(path) as ds:
        if "z" not in ds.variables:
            return None
        z_birth = np.asarray(ds.variables["z"][:], dtype=float).reshape(-1)

    z_birth = z_birth[np.isfinite(z_birth)]
    if z_birth.size == 0:
        return None

    if z_profile.size > 1:
        dz = float(np.median(np.diff(z_profile)))
    else:
        dz = 1.0
    edges = np.concatenate(([z_profile[0] - 0.5 * dz], 0.5 * (z_profile[:-1] + z_profile[1:]), [z_profile[-1] + 0.5 * dz]))
    hist, _ = np.histogram(z_birth, bins=edges)
    return normalize_shape(hist.astype(float))


def mpex_proxy_profiles(args: argparse.Namespace, z_profile: np.ndarray, resonance_z: float) -> dict[str, np.ndarray | str]:
    """Build normalized profile shapes for the MPEX ECH decks.

    If a HERMES axial CSV is supplied, ne/Te/RF-power shapes are interpolated
    from it.  Otherwise a smooth analytic proxy is used.  The output arrays are
    dimensionless shapes; the input deck scalar values set physical units.
    """
    source = "analytic MPEX-like proxy"
    ne_shape = 0.08 + 0.70 * np.exp(-0.5 * np.square((z_profile - 1.7) / 0.55))
    ne_shape += 0.35 * np.exp(-0.5 * np.square((z_profile - 2.8) / 1.75))
    ne_shape += 0.16 / (1.0 + np.exp(-(z_profile - 3.8) / 0.7))
    ne_shape = normalize_shape(ne_shape, args.profile_ne_floor_fraction)

    te_eV = args.profile_te_floor_ev + (args.te_ev - args.profile_te_floor_ev) * (
        0.15 + 0.85 * np.exp(-0.5 * np.square((z_profile - 2.55) / 0.95))
    )
    te_shape = np.clip(te_eV / max(args.te_ev, 1.0e-12), args.profile_te_floor_ev / max(args.te_ev, 1.0e-12), None)
    ti_shape = np.ones_like(z_profile)
    rf_shape = gaussian_shape(z_profile, resonance_z, 0.45)
    source_shape = gaussian_shape(z_profile, args.source_z, args.source_sigma)

    if args.plasma_profile_csv is not None and args.plasma_profile_csv.is_file():
        data = np.genfromtxt(args.plasma_profile_csv, delimiter=",", names=True, dtype=None, encoding=None)
        z_csv = read_csv_column(data, "z_m", "z")
        ne_csv = read_csv_column(data, "Ne_m3", "ne", "density")
        te_csv = read_csv_column(data, "Te_eV", "te", "temperature")
        q_csv = read_csv_column(data, "q_W_m2", "qpar", "q")
        if z_csv is not None:
            order = np.argsort(z_csv)
            z_sorted = z_csv[order]
            if ne_csv is not None:
                ne_interp = np.interp(z_profile, z_sorted, ne_csv[order], left=ne_csv[order][0], right=ne_csv[order][-1])
                ne_shape = normalize_shape(ne_interp, args.profile_ne_floor_fraction)
            if te_csv is not None:
                te_interp = np.interp(z_profile, z_sorted, te_csv[order], left=te_csv[order][0], right=te_csv[order][-1])
                te_interp = np.maximum(te_interp, args.profile_te_floor_ev)
                te_shape = te_interp / max(args.te_ev, 1.0e-12)
            if q_csv is not None:
                q_interp = np.interp(z_profile, z_sorted, q_csv[order], left=0.0, right=q_csv[order][-1])
                rf_shape = normalize_shape(q_interp)
            source = f"HERMES axial CSV {args.plasma_profile_csv.name}"

    if args.source_particles_nc is not None and args.source_particles_nc.is_file():
        particle_source = profile_from_particles_nc(args.source_particles_nc, z_profile)
        if particle_source is not None and np.any(particle_source > 0.0):
            source_shape = particle_source
            source += f"; helicon particle-source z histogram {args.source_particles_nc.name}"
    elif args.use_density_as_source:
        source_shape = normalize_shape(ne_shape)
        source += "; pair source follows density profile"

    return {
        "ne": ne_shape,
        "te": te_shape,
        "ti": ti_shape,
        "source": normalize_shape(source_shape),
        "rf": normalize_shape(rf_shape),
        "description": source,
    }


def fortran_case_path(args: argparse.Namespace) -> Path:
    return args.linear_root / "InputFiles" / f"{args.case_name}.in"


def fortran_output_dir(args: argparse.Namespace) -> Path:
    return args.linear_root / "OutputFiles" / args.case_name / args.descriptor


def picos_input_dir(args: argparse.Namespace) -> Path:
    return args.picos_root / "picosFILES" / "inputFiles"


def picos_hdf_dir(args: argparse.Namespace, tag: str) -> Path:
    return args.picos_root / "picosFILES" / "outputFiles" / tag / "HDF5"


def write_fortran_case(args: argparse.Namespace, z: np.ndarray, b_t: np.ndarray, resonance_z: float) -> Path:
    args.linear_root.mkdir(parents=True, exist_ok=True)
    (args.linear_root / "InputFiles").mkdir(parents=True, exist_ok=True)
    (args.linear_root / "BfieldData").mkdir(parents=True, exist_ok=True)

    bfield_name = f"{args.case_name}_Bfield.txt"
    bfield_path = args.linear_root / "BfieldData" / bfield_name
    np.savetxt(bfield_path, np.column_stack([z, b_t]), fmt="%.16e")

    text = f"""&params_nml
params%fileDescriptor = '{args.descriptor}',
params%repoDir       = "{args.linear_root}",
params%BFieldFileDir = "/BfieldData",
params%BFieldFile    = "/{bfield_name}",
params%nz            = {z.size},

params%NC = {args.particles},
params%NS = {args.steps},
params%dt = {args.physical_time / args.steps:.16e},
params%G  = {args.source_rate:.16e},

params%jstart = {args.steps},
params%jend   = {args.steps},
params%jincr  = 1,

params%dtheta = 6.283185307179586,
params%r1     = 0.0,
params%r2     = {args.radius:.16e},
params%zmax   = {args.z_max:.16e},
params%zmin   = {args.z_min:.16e},
params%NZmesh = {z.size - 2},

params%iSave      = .true.,
params%iPush      = .true.,
params%iColl      = {'.true.' if args.collisions else '.false.'},
params%iHeat      = .true.,
params%iPotential = .false.,
params%iDrag      = .false.,

params%Te0          = {args.te_ev:.16e},
params%Ti0          = {args.ti_ev:.16e},
params%ne0          = {args.ne_m3:.16e},
params%Aion         = 2.0,
params%Zion         = 1.0,
params%species_a    = 1,
params%elevel       = 10,
params%CollOperType = 2,

params%BC_Type    = {args.boundary_type},
params%BC_zp_mean = {args.source_z:.16e},
params%BC_zp_std  = {args.source_sigma:.16e},
params%BC_Ep      = 0.0,
params%BC_Tp      = {args.te_ev:.16e},
params%BC_xip     = 0.707,

params%IC_Type    = {args.ic_type},
params%IC_zp_mean = {args.source_z:.16e},
params%IC_zp_std  = {args.source_sigma:.16e},
params%IC_Ep      = 0.0,
params%IC_Tp      = {args.te_ev:.16e},
params%IC_xip     = 0.707,

params%f_RF       = {args.rf_frequency:.16e},
params%zRes1      = {resonance_z - args.rf_window_half_width:.16e},
params%zRes2      = {args.rf_window_end:.16e},
params%kper       = {args.kper:.16e},
params%kpar       = {args.kpar:.16e},
params%Prf        = {args.rf_power:.16e},
params%n_harmonic = {args.n_harmonic},

params%s1   = 0.0,
params%s2   = 1.3,
params%s3   = 4.0,
params%phi1 = 0.0,
params%phi2 = 0.0,
params%phi3 = 0.0
/
"""
    path = fortran_case_path(args)
    path.write_text(text)
    return path


def rf_block(args: argparse.Namespace, rf_file: str, simulation_time_gyro: float, resonance_z: float) -> str:
    return f"""// Ion RF operator:
// =============================================================================
RF_ion_Prf                      0.0
RF_ion_n_harmonic               {args.n_harmonic}
RF_ion_freq                     {args.rf_frequency:.16e}
RF_ion_x1                       {resonance_z - args.rf_window_half_width:.16e}
RF_ion_x2                       {args.rf_window_end:.16e}
RF_ion_t_ON                     0.0
RF_ion_t_OFF                    {simulation_time_gyro:.16e}
RF_ion_kpar                     {args.kpar:.16e}
RF_ion_kper                     {args.kper:.16e}
RF_ion_handedness               -1
RF_ion_EfieldMode               0
RF_ion_resonanceMode            {args.rf_resonance_mode}
RF_ion_EfieldAmplitude          {args.rf_efield_amplitude:.16e}
RF_ion_maxEnergyGainFraction    0.0
RF_ion_maxParticleEnergy        0.0
RF_ion_maxVelocityFractionC     0.0
RF_ion_Prf_fileName             {rf_file}
RF_ion_Prf_NS                   {args.profile_points}

// Electron RF/ECH operator:
// =============================================================================
RF_electron_Prf                      {args.rf_power:.16e}
RF_electron_n_harmonic               {args.n_harmonic}
RF_electron_freq                     {args.rf_frequency:.16e}
RF_electron_x1                       {resonance_z - args.rf_window_half_width:.16e}
RF_electron_x2                       {args.rf_window_end:.16e}
RF_electron_t_ON                     0.0
RF_electron_t_OFF                    {simulation_time_gyro:.16e}
RF_electron_kpar                     {args.kpar:.16e}
RF_electron_kper                     {args.kper:.16e}
RF_electron_handedness               -1
RF_electron_EfieldMode               0
RF_electron_resonanceMode            {args.rf_resonance_mode}
RF_electron_EfieldAmplitude          {args.rf_efield_amplitude:.16e}
RF_electron_maxEnergyGainFraction    {args.rf_max_energy_gain_fraction:.16e}
RF_electron_maxParticleEnergy        {args.rf_max_particle_energy:.16e}
RF_electron_maxVelocityFractionC     {args.rf_max_velocity_fraction_c:.16e}
RF_electron_Prf_fileName             {rf_file}
RF_electron_Prf_NS                   {args.profile_points}

"""


def write_picos_case(args: argparse.Namespace, tag: str, z_b: np.ndarray, b_norm: np.ndarray, resonance_z: float, relativistic: bool) -> tuple[Path, Path]:
    input_dir = picos_input_dir(args)
    input_dir.mkdir(parents=True, exist_ok=True)
    b_norm_file = f"{tag}_B_norm.txt"
    ne_file = f"{tag}_ne_norm.txt"
    te_file = f"{tag}_te_norm.txt"
    ti_file = f"{tag}_ti_norm.txt"
    source_file = f"{tag}_pair_source_norm.txt"
    rf_file = f"{tag}_rf_power_norm.txt"
    z_b_profile = picos_profile_z(args)
    z_aux_profile = picos_aux_profile_z(args)
    b_for_picos = np.interp(z_b_profile, z_b, b_norm, left=b_norm[0], right=b_norm[-1])
    profiles = mpex_proxy_profiles(args, z_aux_profile, resonance_z)
    write_vector(input_dir / b_norm_file, b_for_picos)
    write_vector(input_dir / ne_file, profiles["ne"])
    write_vector(input_dir / te_file, profiles["te"])
    write_vector(input_dir / ti_file, profiles["ti"])
    write_vector(input_dir / source_file, profiles["source"])
    write_vector(input_dir / rf_file, profiles["rf"])

    cv_b = args.cv_b
    dx = (args.z_max - args.z_min) / (args.profile_points - 2)
    dp = dx / ion_skin_depth(args.ne_m3)
    simulation_time_gyro = args.physical_time / ion_gyroperiod(cv_b)
    output_count = max(1, int(args.output_count))
    output_cadence_gyro = simulation_time_gyro / output_count
    electron_ppc = max(1, int(round(args.particles / max(args.profile_points - 2, 1))))
    ion_ppc = max(1, int(round(electron_ppc / 4.0)))

    field_description = (
        "reformulated Poisson electrostatic field solve"
        if args.efield_solve and args.field_solve_model == 2
        else "standard electrostatic Poisson field solve"
        if args.efield_solve and args.field_solve_model == 1
        else "electrostatic solve disabled for direct Fortran RF-transport comparison"
    )

    input_text = f"""// PICOS input generated for MPEX scenario-14 Fig. 17 E-z comparison
// Magnetic field source profile: templateFILES/MPEX_B_norm_PICOS_scenario_14.txt.
// Plasma/source profile model: {profiles["description"]}.
// Model: 1D-2V guiding-center D+ ions plus kinetic electrons.
// Heating: electron ECH enabled with scenario-14 B-field; ion RF heating disabled.
// Field: {field_description}.
// Profile files are dimensionless shapes; IC_ne, IC_Te, IC_Tpar, and IC_Tper set physical amplitudes.
// RF_*_Prf_fileName is written for bookkeeping/future GENRAY coupling; the current RF operator uses scalar RF_*_Prf in power-balance mode.
// =============================================================================
mpisForFields               {args.mpis_for_fields}
quietStart                  {args.quiet_start}
IC_velocityDistributionModel 1
IC_randomSeed              271828
IC_weightScale              {args.ic_weight_scale:.16e}
numberOfParticleSpecies     2
numberOfTracerSpecies       0
advanceParticleMethod       1

// Characteristic values:
// =============================================================================
CV_ne                       {args.ne_m3:.16e}
CV_Te                       {args.te_ev:.16e}
CV_B                        {cv_b:.16e}
CV_Tpar                     {args.ti_ev:.16e}
CV_Tper                     {args.ti_ev:.16e}

// Simulation time:
// =============================================================================
DTc                         {args.dtc:.16e}
simulationTime              {simulation_time_gyro:.16e}

// Switches:
// =============================================================================
SW_EfieldSolve              {args.efield_solve}
SW_fieldSolveModel          {args.field_solve_model}
SW_electronGyroTimeStepLimiter   {args.electron_gyro_timestep_limiter}
SW_electronPlasmaTimeStepLimiter {args.electron_plasma_timestep_limiter}
SW_BfieldSolve              0
SW_Collisions               {1 if args.collisions else 0}
CollOperType                2
SW_collisionConservationProjection 0
collisionRandomSeed        314159
SW_RFheating                {args.rf_heating}
SW_RFheatingIons            0
SW_RFheatingElectrons       1
SW_relativisticElectrons    {1 if relativistic else 0}
SW_relativisticRFOperator   {1 if relativistic else 0}
SW_pairSource               {args.pair_source}
SW_advancePos               1
SW_linearSolve              0

// Restart controls:
// restart_snapshot=-1 loads the latest numeric HDF5 output snapshot.
// For a two-stage run, first run SW_RFheating=0 to steady state, then set
// restart_enabled=1 and restart_path to the stage-1 output HDF5 directory.
// =============================================================================
restart_enabled             {1 if args.restart_path is not None else 0}
restart_path                {args.restart_path if args.restart_path is not None else "none"}
restart_snapshot            {args.restart_snapshot}
restart_continueTime        {1 if args.restart_continue_time else 0}
restart_particleFilePrefix  PARTICLES_FILE_
restart_fieldsFilePrefix    FIELDS_FILE_

// Magnetic field initial conditions:
// =============================================================================
IC_uniformBfield            0
IC_BX                       {cv_b:.16e}
IC_BY                       0.0
IC_BZ                       0.0
IC_BX_NX                    {args.profile_points}
IC_BX_fileName              {b_norm_file}
IC_phiLeft                  0.0
IC_phiRight                 0.0
Poisson_BCModel             1
Poisson_sheathCoefficient   3.0
ReformulatedPoisson_lambda  -1.0
ReformulatedPoisson_quasiNeutral {args.reformulated_poisson_quasineutral}

// Geometry:
// =============================================================================
dp                          {dp:.16e}
r1                          0.0
r2                          {args.radius:.16e}
LX_min                      {args.z_min:.16e}
LX_max                      {args.z_max:.16e}

// Electron initial conditions:
// =============================================================================
IC_ne                       {args.ne_m3:.16e}
IC_Te                       {args.te_ev:.16e}
IC_Te_NX                    {args.profile_points}
IC_Te_fileName              {te_file}

// Coupled electron-ion source:
// =============================================================================
pairSource_ionSpecies       1
pairSource_electronSpecies  2
pairSource_rate             {args.source_rate:.16e}
pairSource_mean_x           {args.source_z:.16e}
pairSource_sigma_x          {args.source_sigma:.16e}
pairSource_Ti_birth         {args.ti_ev:.16e}
pairSource_Te_birth         {args.te_ev:.16e}
pairSource_Ei_birth         0
pairSource_Ee_birth         0
pairSource_eta_i            0
pairSource_eta_e            0
pairSource_positionMode     0
pairSource_weightMode       {args.pair_source_weight_mode}
pairSource_fileName         {source_file}
pairSource_NS               {args.profile_points}
pairSource_maxParticleWeight 1000

{rf_block(args, rf_file, simulation_time_gyro, resonance_z)}
// Output variables:
// =============================================================================
// outputCadence is in background ion gyroperiod units. This deck requests
// {output_count} output intervals plus the initial t=0 snapshot.
outputCadence               {output_cadence_gyro:.16e}
outputs_variables           {{X_p,V_p,a_p,mu_p,BX_p,EX_p,BX_m,dBX_m,ddBX_m,n_m,Tpar_m,Tper_m,u_m,EX_m,Phi_m}}

// Data smoothing:
// =============================================================================
smoothingParameter          1.0000000000000000e-03
filtersPerIterationFields   2
filtersPerIterationIons     2
"""

    ions_text = f"""// PICOS species generated for MPEX scenario-14 Fig. 17 E-z comparison
// Species are 1D-2V guiding-center D+ ions and kinetic electrons.
// =============================================================================
// Species 1: deuterium ions
// =============================================================================
SPECIES1                      1
NPC1                          {ion_ppc}
pctSupPartOutput1             1.0000000000000000e+02
Z1                            1
M1                            2.0000000000000000e+00

IC_type_1                     1
IC_Tper_1                     {args.ti_ev:.16e}
IC_Tper_fileName_1            {ti_file}
IC_Tper_NX_1                  {args.profile_points}
IC_Tpar_1                     {args.ti_ev:.16e}
IC_Tpar_fileName_1            {ti_file}
IC_Tpar_NX_1                  {args.profile_points}
IC_densityFraction_1          1.0
IC_densityFraction_fileName_1 {ne_file}
IC_densityFraction_NX_1       {args.profile_points}
IC_weightScale_1              {args.ic_weight_scale:.16e}

BC_type_1                     {args.boundary_type}
BC_T_1                        {args.ti_ev:.16e}
BC_E_1                        0.0
BC_eta_1                      0.7853981633974483
BC_mean_x_1                   {args.source_z:.16e}
BC_sigma_x_1                  {args.source_sigma:.16e}
BC_G_1                        {args.source_rate:.16e}
BC_G_fileName_1               {source_file}
BC_G_NS_1                     {args.profile_points}

// =============================================================================
// Species 2: kinetic electrons
// =============================================================================
SPECIES2                      1
NPC2                          {electron_ppc}
pctSupPartOutput2             1.0000000000000000e+02
Z2                            -1
M2                            5.485799090441e-4

IC_type_2                     1
IC_Tper_2                     {args.te_ev:.16e}
IC_Tper_fileName_2            {te_file}
IC_Tper_NX_2                  {args.profile_points}
IC_Tpar_2                     {args.te_ev:.16e}
IC_Tpar_fileName_2            {te_file}
IC_Tpar_NX_2                  {args.profile_points}
IC_densityFraction_2          1.0
IC_densityFraction_fileName_2 {ne_file}
IC_densityFraction_NX_2       {args.profile_points}
IC_weightScale_2              {args.ic_weight_scale:.16e}

BC_type_2                     {args.boundary_type}
BC_T_2                        {args.te_ev:.16e}
BC_E_2                        0.0
BC_eta_2                      0.7853981633974483
BC_mean_x_2                   {args.source_z:.16e}
BC_sigma_x_2                  {args.source_sigma:.16e}
BC_G_2                        {args.source_rate:.16e}
BC_G_fileName_2               {source_file}
BC_G_NS_2                     {args.profile_points}
"""
    input_path = input_dir / f"input_file_{tag}.input"
    ions_path = input_dir / f"ions_properties_{tag}.ion"
    input_path.write_text(input_text)
    ions_path.write_text(ions_text)
    return input_path, ions_path


def setup_cases(args: argparse.Namespace) -> None:
    z, b_t, b_norm = scenario14_profile(args)
    args.profile_points = z.size
    args.cv_b = float(np.max(np.abs(b_t)))
    resonance_z = selected_resonance_z(args, z, b_t)
    args.rf_window_end = args.rf_window_end if args.rf_window_end is not None else args.z_max
    write_fortran_case(args, z, b_t, resonance_z)
    if args.picos_run_mode in ("both", "nonrel"):
        write_picos_case(args, args.picos_tag_nonrel, z, b_norm, resonance_z, relativistic=False)
    if args.picos_run_mode in ("both", "rel"):
        write_picos_case(args, args.picos_tag_rel, z, b_norm, resonance_z, relativistic=True)


def read_fortran_snapshot(args: argparse.Namespace) -> dict[str, np.ndarray] | None:
    output_dir = fortran_output_dir(args)
    if not (output_dir / "zp.out").is_file() or not (output_dir / "kep.out").is_file():
        return None
    z = read_fortran_record_float64(output_dir / "zp.out")
    energy = read_fortran_record_float64(output_dir / "kep.out")
    cols = max(1, z.size // args.particles)
    z = z.reshape((args.particles, cols), order="F")[:, -1]
    energy = energy.reshape((args.particles, cols), order="F")[:, -1]
    return {"z_m": z, "energy_eV": energy}


def kinetic_energy_eV(speed: np.ndarray, relativistic: bool) -> np.ndarray:
    if not relativistic:
        return 0.5 * M_E * speed * speed / E_CHARGE
    beta2 = np.square(speed / C_LIGHT)
    gamma = 1.0 / np.sqrt(1.0 - np.clip(beta2, 0.0, 1.0 - 1.0e-15))
    return (gamma - 1.0) * M_E * C_LIGHT * C_LIGHT / E_CHARGE


def read_picos_snapshot(args: argparse.Namespace, tag: str, relativistic: bool) -> dict[str, Any] | None:
    hdf_dir = picos_hdf_dir(args, tag)
    files = sorted(hdf_dir.glob("PARTICLES_FILE_*.h5"))
    if not files:
        return None
    positions = []
    velocities = []
    rf_power = math.nan
    for path in files:
        try:
            components = picos_components(read_hdf5_dataset(path, "/1/ions/spp_2/V_p"))
            x = read_hdf5_dataset(path, "/1/ions/spp_2/X_p").reshape(-1)
        except Exception:
            continue
        velocities.append(components)
        positions.append(x)
        if not math.isfinite(rf_power):
            rf_power = read_hdf5_scalar_optional(path, "/1/rf/electron/E3")
    if not velocities:
        return None
    v = np.concatenate(velocities, axis=1)
    speed = np.sqrt(np.sum(v * v, axis=0))
    return {"z_m": np.concatenate(positions), "energy_eV": kinetic_energy_eV(speed, relativistic), "rf_absorbed_power_W": rf_power}


def gaussian_kernel1d(sigma_bins: float) -> np.ndarray:
    if sigma_bins <= 0.0:
        return np.ones(1)
    radius = max(1, int(math.ceil(3.0 * sigma_bins)))
    x = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-0.5 * np.square(x / sigma_bins))
    return kernel / kernel.sum()


def convolve_axis_reflect(values: np.ndarray, kernel: np.ndarray, axis: int) -> np.ndarray:
    if kernel.size == 1:
        return values
    radius = kernel.size // 2
    moved = np.moveaxis(values, axis, 0)
    padded = np.pad(moved, [(radius, radius), *[(0, 0)] * (moved.ndim - 1)], mode="reflect")
    out = np.empty_like(moved, dtype=float)
    for idx in np.ndindex(moved.shape[1:]):
        out[(slice(None),) + idx] = np.convolve(padded[(slice(None),) + idx], kernel, mode="valid")
    return np.moveaxis(out, 0, axis)


def histogram_log_density(z: np.ndarray, energy: np.ndarray, args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mask = np.isfinite(z) & np.isfinite(energy) & (energy >= 0.0)
    hist, z_edges, e_edges = np.histogram2d(
        z[mask],
        energy[mask],
        bins=(args.z_bins, args.energy_bins),
        range=((args.plot_z_min, args.plot_z_max), (args.energy_min, args.energy_max)),
    )
    hist = convolve_axis_reflect(hist, gaussian_kernel1d(args.smooth_z_bins), axis=0)
    hist = convolve_axis_reflect(hist, gaussian_kernel1d(args.smooth_energy_bins), axis=1)
    if hist.max() > 0.0:
        density = hist / hist.sum()
        log_density = np.log10(np.where(density > 0.0, density, np.nan))
        if args.mask_below_log_vmin:
            log_density = np.where(log_density >= args.log_vmin, log_density, np.nan)
    else:
        log_density = np.full_like(hist, np.nan, dtype=float)
    return z_edges, e_edges, log_density.T


def plot_single_map(
    args: argparse.Namespace,
    label: str,
    data: dict[str, Any],
    z_b: np.ndarray,
    b_t: np.ndarray,
    resonance_z: float,
    output_path: Path,
) -> None:
    plt = ensure_plotting()
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    z_edges, e_edges, log_density = histogram_log_density(data["z_m"], data["energy_eV"], args)
    fig, ax = plt.subplots(figsize=(13.2, 7.2), constrained_layout=True)
    mesh = ax.pcolormesh(
        z_edges,
        e_edges,
        log_density,
        cmap=args.colormap,
        shading="auto",
        vmin=args.log_vmin,
        vmax=args.log_vmax,
    )
    ax.set_xlim(args.plot_z_min, args.plot_z_max)
    ax.set_ylim(args.energy_min, args.energy_max)
    ax.set_xlabel("z [m]")
    ax.set_ylabel("E_e [eV]")
    ax.grid(True, alpha=0.25)
    title = args.plot_title if not label else f"{args.plot_title} ({label})"
    ax.set_title(title)
    ax.axvline(resonance_z, color="red", ls=":", lw=2.5)
    ax.axvline(args.target_z, color="limegreen", lw=3.0)

    def b_to_energy_axis(b_value: np.ndarray | float) -> np.ndarray | float:
        return args.energy_min + (np.asarray(b_value) - args.b_axis_min) * (
            (args.energy_max - args.energy_min) / (args.b_axis_max - args.b_axis_min)
        )

    def energy_axis_to_b(energy_value: np.ndarray | float) -> np.ndarray | float:
        return args.b_axis_min + (np.asarray(energy_value) - args.energy_min) * (
            (args.b_axis_max - args.b_axis_min) / (args.energy_max - args.energy_min)
        )

    ax.plot(z_b, b_to_energy_axis(b_t), color="red", lw=2.2)
    ax2 = ax.secondary_yaxis("right", functions=(energy_axis_to_b, b_to_energy_axis))
    ax2.set_ylabel("B0 [T]", color="red")
    ax2.tick_params(axis="y", colors="red")
    cax = inset_axes(ax, width="2.5%", height="74%", loc="upper right", borderpad=2.1)
    cbar = fig.colorbar(mesh, cax=cax)
    cbar.ax.set_title("log10(f(E,z))", fontsize=27, pad=6)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_comparison(args: argparse.Namespace) -> list[Path]:
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    z, b_t, _ = scenario14_profile(args)
    resonance_z = selected_resonance_z(args, z, b_t)
    args.rf_window_end = args.rf_window_end if args.rf_window_end is not None else args.z_max

    outputs = {
        "fortran_nonrel": read_fortran_snapshot(args),
        "picos_nonrel": read_picos_snapshot(args, args.picos_tag_nonrel, relativistic=False),
        "picos_rel": read_picos_snapshot(args, args.picos_tag_rel, relativistic=True),
    }
    labels = {
        "fortran_nonrel": "Fortran",
        "picos_nonrel": "PICOS++ nonrel",
        "picos_rel": "PICOS++ rel",
    }
    paths: list[Path] = []
    for name, data in outputs.items():
        if data is None:
            continue
        path = out_dir / f"{name}_fig17_Ez_map.png"
        plot_single_map(args, labels[name], data, z, b_t, resonance_z, path)
        paths.append(path)

    plot_b_profile(args, z, b_t, resonance_z)
    return paths


def plot_b_profile(args: argparse.Namespace, z: np.ndarray, b_t: np.ndarray, resonance_z: float) -> Path:
    plt = ensure_plotting()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    path = args.out_dir / "scenario14_bfield_overlay.png"
    fig, ax = plt.subplots(figsize=(13.2, 4.8), constrained_layout=True)
    ax.plot(z, b_t, color="red", lw=2.3, label="scenario-14 B0")
    ax.axhline(resonance_b(args), color="0.25", ls="--", lw=1.2, label=f"{args.n_harmonic}nd harmonic B")
    ax.axvline(resonance_z, color="red", ls=":", lw=2.0, label="selected resonance")
    ax.axvline(args.target_z, color="limegreen", lw=2.5, label="target")
    if abs(args.rf_window_end - args.target_z) > 1.0e-12:
        ax.axvline(args.rf_window_end, color="orange", ls="--", lw=1.5, label="RF window end")
    ax.set_xlim(args.plot_z_min, args.plot_z_max)
    ax.set_ylim(args.b_axis_min, args.b_axis_max)
    ax.set_xlabel("z [m]")
    ax.set_ylabel("B0 [T]")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=24)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def write_summary(args: argparse.Namespace, plot_paths: list[Path]) -> Path:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    z, b_t, _ = scenario14_profile(args)
    resonance_z = selected_resonance_z(args, z, b_t)
    rows = []
    for name, data in (
        ("fortran_nonrel", read_fortran_snapshot(args)),
        ("picos_nonrel", read_picos_snapshot(args, args.picos_tag_nonrel, relativistic=False)),
        ("picos_rel", read_picos_snapshot(args, args.picos_tag_rel, relativistic=True)),
    ):
        if data is None:
            rows.append({"run": name, "has_output": False})
            continue
        energy = data["energy_eV"]
        rows.append(
            {
                "run": name,
                "has_output": True,
                "count": int(np.isfinite(energy).sum()),
                "mean_energy_eV": float(np.nanmean(energy)),
                "p99_energy_eV": float(np.nanpercentile(energy, 99)),
                "max_energy_eV": float(np.nanmax(energy)),
                "rf_absorbed_power_W": float(data.get("rf_absorbed_power_W", math.nan)),
            }
        )

    csv_path = args.out_dir / "scenario14_fig17_summary.csv"
    fieldnames = ["run", "has_output", "count", "mean_energy_eV", "p99_energy_eV", "max_energy_eV", "rf_absorbed_power_W"]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    path = args.out_dir / "README.md"
    lines = [
        "# MPEX Scenario-14 Fig. 17 E-z Comparison",
        "",
        "This workflow targets the figure shown in the screenshot: electron kinetic energy versus axial coordinate, colored by `log10(f(E,z))`, with the MPEX scenario-14 magnetic-field profile overlaid.",
        "",
        "## Inputs",
        "",
        f"- Scenario-14 B template: `{args.scenario14_b_file}`.",
        f"- Fortran input generated: `{fortran_case_path(args)}`.",
        f"- PICOS++ nonrel input generated: `{picos_input_dir(args) / ('input_file_' + args.picos_tag_nonrel + '.input')}`.",
        f"- PICOS++ relativistic input generated: `{picos_input_dir(args) / ('input_file_' + args.picos_tag_rel + '.input')}`.",
        "",
        "## Key Settings",
        "",
        f"- Domain used for scenario-14 active profile: `{args.z_min:g}` to `{args.z_max:g}` m.",
        f"- Plot axis: `{args.plot_z_min:g}` to `{args.plot_z_max:g}` m and `{args.energy_min:g}` to `{args.energy_max:g}` eV.",
        f"- B scale applied before normalization: `{args.b_scale:g}` T, giving `Bmax={float(np.max(b_t)):.4g}` T.",
        f"- {args.n_harmonic}nd harmonic resonance field for `{args.rf_frequency:.4g}` Hz: `{resonance_b(args):.4g}` T.",
        f"- Maximum allowed RF power in this comparison: `{args.max_rf_power:.4g}` W.",
        f"- Selected resonance marker: `z={resonance_z:.4g}` m.",
        f"- Target marker: `z={args.target_z:.4g}` m.",
        f"- RF window: `z={resonance_z - args.rf_window_half_width:.4g}` m to `z={args.rf_window_end:.4g}` m.",
        f"- RF resonance mode: `{args.rf_resonance_mode}`.",
        f"- RF heating switch: `{args.rf_heating}`.",
        f"- Restart path: `{args.restart_path}`.",
        f"- Restart snapshot: `{args.restart_snapshot}`.",
        f"- Restart continues physical clock: `{args.restart_continue_time}`.",
        f"- Fortran/PICOS physical time requested: `{args.physical_time:.4g}` s.",
        f"- PICOS output intervals requested: `{args.output_count}` plus the initial snapshot.",
        "",
        "## Outputs",
        "",
        f"- Summary table: `{csv_path.name}`.",
        f"- B overlay: `scenario14_bfield_overlay.png`.",
    ]
    for plot in plot_paths:
        lines.append(f"- E-z map: `{plot.name}`.")
    lines.extend(
        [
            "",
            "## Current Status",
            "",
            "| run | output | count | mean E [eV] | P99 [eV] | max E [eV] | RF absorbed [W] |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in rows:
        if not row["has_output"]:
            lines.append(f"| {row['run']} | no |  |  |  |  |  |")
            continue
        lines.append(
            f"| {row['run']} | yes | {row['count']} | {row['mean_energy_eV']:.4g} | "
            f"{row['p99_energy_eV']:.4g} | {row['max_energy_eV']:.4g} | {row['rf_absorbed_power_W']:.4g} |"
        )
    lines.extend(
        [
            "",
            "## Commands",
            "",
            "```bash",
            "cd /Users/78k/Desktop/picos_kinetic_electron_eval/PICOS",
            "PYTHONPYCACHEPREFIX=/private/tmp/picos_pycache \\",
            "python3 scripts/compare_mpex_scenario14_fig17.py --setup --plot",
            "",
            "PYTHONPYCACHEPREFIX=/private/tmp/picos_pycache MPLCONFIGDIR=/private/tmp/picos_mpl \\",
            "python3 scripts/compare_mpex_scenario14_fig17.py --run-fortran --plot",
            "",
            "# Full PICOS++ physical-time run is intentionally opt-in and may be very expensive.",
            "PYTHONPYCACHEPREFIX=/private/tmp/picos_pycache MPLCONFIGDIR=/private/tmp/picos_mpl \\",
            "python3 scripts/compare_mpex_scenario14_fig17.py --run-picos --plot",
            "```",
        ]
    )
    path.write_text("\n".join(lines) + "\n")
    return path


def run_fortran(args: argparse.Namespace) -> None:
    import os

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(args.fortran_omp_threads or args.omp_threads)
    env["OMP_PROC_BIND"] = "false"
    env["REPO_DIR"] = str(args.linear_root)
    env["INPUT_FILE"] = f"{args.case_name}.in"
    env["INPUT_FILE_DIR"] = str(fortran_case_path(args))
    subprocess.check_call([str(args.linear_root / "src" / "linFP")], cwd=args.linear_root, env=env)


def run_picos(args: argparse.Namespace) -> None:
    import os

    env = os.environ.copy()
    env["HDF5_USE_FILE_LOCKING"] = "FALSE"
    env["OMP_NUM_THREADS"] = str(args.picos_omp_threads or args.omp_threads)
    xpicos = args.picos_build_dir / "picosFILES" / "src" / "xpicos"
    picos_files = args.picos_root / "picosFILES"
    tags = []
    if args.picos_run_mode in ("both", "nonrel"):
        tags.append(args.picos_tag_nonrel)
    if args.picos_run_mode in ("both", "rel"):
        tags.append(args.picos_tag_rel)
    for tag in tags:
        command = [str(xpicos), "1-D", "outputFiles", tag]
        if args.mpi_ranks > 1:
            command = ["mpirun", "-np", str(args.mpi_ranks), *command]
        subprocess.check_call(command, cwd=picos_files, env=env)


def normalize_args(args: argparse.Namespace) -> argparse.Namespace:
    args.picos_root = args.picos_root.resolve()
    args.linear_root = args.linear_root.resolve()
    if not args.picos_build_dir.is_absolute():
        args.picos_build_dir = args.picos_root / args.picos_build_dir
    args.picos_build_dir = args.picos_build_dir.resolve()
    args.scenario14_b_file = args.scenario14_b_file.resolve()
    if args.plasma_profile_csv is not None:
        args.plasma_profile_csv = args.plasma_profile_csv.expanduser().resolve()
    if args.source_particles_nc is not None:
        args.source_particles_nc = args.source_particles_nc.expanduser().resolve()
    if args.restart_path is not None:
        args.restart_path = args.restart_path.expanduser()
        if not args.restart_path.is_absolute():
            args.restart_path = args.restart_path.resolve()
    if not args.out_dir.is_absolute():
        args.out_dir = args.picos_root / args.out_dir
    args.out_dir = args.out_dir.resolve()
    if args.rf_window_end is None:
        args.rf_window_end = args.z_max
    if args.target_z is None:
        args.target_z = args.z_max
    raw_b = np.loadtxt(args.scenario14_b_file)
    args.profile_points = int(raw_b.size)
    if args.b_scale is None:
        if args.resonance_b_t is None:
            raise ValueError("--b-scale is required when --resonance-b-t is not set")
        if args.scale_reference_z_to_resonance_b:
            z_raw = np.linspace(args.z_min, args.z_max, raw_b.size)
            b_ref = float(np.interp(args.reference_resonance_z, z_raw, raw_b))
            if b_ref <= 0.0:
                raise ValueError("--reference-resonance-z maps to a non-positive B template value")
            args.b_scale = float(args.resonance_b_t) / b_ref
        else:
            args.b_scale = float(args.resonance_b_t) / float(np.max(np.abs(raw_b)))
    if args.resonance_b_t is not None and args.derive_rf_frequency_from_resonance_b:
        args.rf_frequency = rf_frequency_for_resonance_b(float(args.resonance_b_t), args.n_harmonic)
    if args.rf_power > args.max_rf_power * (1.0 + 1.0e-12):
        raise ValueError(f"RF power {args.rf_power:g} W exceeds the configured cap {args.max_rf_power:g} W")
    if args.mpi_ranks < 2 or args.mpi_ranks % 2 != 0:
        raise ValueError("PICOS++ requires an even number of MPI ranks; use --mpi-ranks 2 for the quick comparison")
    _, b_t, _ = scenario14_profile(args)
    args.cv_b = float(np.max(np.abs(b_t)))
    if args.b_axis_max is None:
        args.b_axis_max = max(float(np.max(b_t)) * 1.05, resonance_b(args) * 1.05)
    return args


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--picos-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--picos-build-dir", type=Path, default=Path("build"), help="PICOS++ build directory, absolute or relative to --picos-root.")
    parser.add_argument("--linear-root", type=Path, default=Path("/Users/78k/Desktop/picos_kinetic_electron_eval/LinearFokkerPlanck_Axisymmetric"))
    parser.add_argument("--scenario14-b-file", type=Path, default=Path("templateFILES/MPEX_B_norm_PICOS_scenario_14.txt"))
    parser.add_argument("--out-dir", type=Path, default=Path("validation/mpex_scenario14_fig17"))
    parser.add_argument("--case-name", default="xp_MPEX_Scenario14_Fig17_B120_Target8_300kW")
    parser.add_argument("--descriptor", default="scenario14_fig17_b120_target8_300kW")
    parser.add_argument("--picos-tag-nonrel", default="mpex_scenario14_fig17_b120_target8_300kW_nonrel")
    parser.add_argument("--picos-tag-rel", default="mpex_scenario14_fig17_b120_target8_300kW_rel")
    parser.add_argument("--particles", type=int, default=2048)
    parser.add_argument("--steps", type=int, default=10000)
    parser.add_argument("--physical-time", type=float, default=2.0e-5)
    parser.add_argument("--output-count", type=int, default=1, help="Number of PICOS output intervals over the requested physical time; the initial t=0 snapshot is also written.")
    parser.add_argument("--ne-m3", type=float, default=1.0e19)
    parser.add_argument("--te-ev", type=float, default=15.0)
    parser.add_argument("--ti-ev", type=float, default=15.0)
    parser.add_argument("--b-scale", type=float, default=None, help="Multiplier applied to the checked-in scenario-14 B profile before normalization.")
    parser.add_argument("--scale-reference-z-to-resonance-b", action=argparse.BooleanOptionalAction, default=True, help="When --b-scale is omitted, scale B so B(reference resonance z) equals --resonance-b-t. Disable to scale the profile peak instead.")
    parser.add_argument("--z-min", type=float, default=-2.0)
    parser.add_argument("--z-max", type=float, default=8.0)
    parser.add_argument("--plot-z-min", type=float, default=-2.0)
    parser.add_argument("--plot-z-max", type=float, default=8.0)
    parser.add_argument("--plot-title", default="I: 1900 [A]")
    parser.add_argument("--energy-min", type=float, default=0.0)
    parser.add_argument("--energy-max", type=float, default=3500.0)
    parser.add_argument("--z-bins", type=int, default=160)
    parser.add_argument("--energy-bins", type=int, default=160)
    parser.add_argument("--smooth-z-bins", type=float, default=0.0, help="Visualization-only Gaussian smoothing width in z histogram bins.")
    parser.add_argument("--smooth-energy-bins", type=float, default=0.0, help="Visualization-only Gaussian smoothing width in energy histogram bins.")
    parser.add_argument("--log-vmin", type=float, default=-6.0)
    parser.add_argument("--log-vmax", type=float, default=-2.0)
    parser.add_argument("--mask-below-log-vmin", action=argparse.BooleanOptionalAction, default=False, help="Hide histogram cells below --log-vmin instead of drawing them with the under-range colormap.")
    parser.add_argument("--b-axis-min", type=float, default=0.0)
    parser.add_argument("--b-axis-max", type=float, default=None)
    parser.add_argument("--colormap", default="viridis")
    parser.add_argument("--rf-frequency", type=float, default=28.0e9)
    parser.add_argument("--n-harmonic", type=int, default=2)
    parser.add_argument("--resonance-b-t", type=float, default=1.2, help="Target ECH resonance magnetic field in tesla. By default RF frequency is derived from this value.")
    parser.add_argument("--derive-rf-frequency-from-resonance-b", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--rf-power", type=float, default=3.0e5)
    parser.add_argument("--rf-heating", type=int, choices=[0, 1], default=1)
    parser.add_argument("--max-rf-power", type=float, default=3.0e5)
    parser.add_argument("--rf-efield-amplitude", type=float, default=1.0e4)
    parser.add_argument("--rf-max-energy-gain-fraction", type=float, default=0.25)
    parser.add_argument("--rf-max-particle-energy", type=float, default=5000.0)
    parser.add_argument("--rf-max-velocity-fraction-c", type=float, default=0.2)
    parser.add_argument("--rf-resonance-mode", type=int, choices=[0, 1], default=0, help="0 uses sign crossing; 1 uses Fortran-compatible resNum<0 window flag.")
    parser.add_argument("--efield-solve", type=int, choices=[0, 1], default=0)
    parser.add_argument("--field-solve-model", type=int, choices=[0, 1, 2], default=0)
    parser.add_argument("--electron-gyro-timestep-limiter", type=int, choices=[0, 1], default=0)
    parser.add_argument("--electron-plasma-timestep-limiter", type=int, choices=[0, 1], default=0)
    parser.add_argument("--reformulated-poisson-quasineutral", type=int, choices=[0, 1], default=0)
    parser.add_argument("--kpar", type=float, default=5.864e3)
    parser.add_argument("--kper", type=float, default=2.992e4)
    parser.add_argument("--reference-resonance-z", type=float, default=3.1)
    parser.add_argument("--resonance-z", type=float, default=None)
    parser.add_argument("--rf-window-half-width", type=float, default=0.02)
    parser.add_argument("--rf-window-end", type=float, default=None, help="End of active RF heating window in meters; defaults to --z-max.")
    parser.add_argument("--target-z", type=float, default=None, help="Target marker position in meters; defaults to --z-max.")
    parser.add_argument("--boundary-type", type=int, default=1)
    parser.add_argument("--source-rate", type=float, default=1.0e23)
    parser.add_argument("--source-z", type=float, default=0.0)
    parser.add_argument("--source-sigma", type=float, default=0.15)
    parser.add_argument("--plasma-profile-csv", type=Path, default=None, help="Optional axial profile CSV with z_m, Ne_m3, Te_eV, and optionally q_W_m2 columns.")
    parser.add_argument("--source-particles-nc", type=Path, default=None, help="Optional helicon source-particle NetCDF; the z histogram is used as the pair-source profile.")
    parser.add_argument("--profile-ne-floor-fraction", type=float, default=0.05, help="Minimum density-shape value as a fraction of the profile peak.")
    parser.add_argument("--profile-te-floor-ev", type=float, default=0.5, help="Minimum Te used when building normalized Te profile files.")
    parser.add_argument("--use-density-as-source", action=argparse.BooleanOptionalAction, default=False, help="Use the density profile shape as the pair-source shape when no source particle file is provided.")
    parser.add_argument("--quiet-start", type=int, choices=[0, 1], default=0)
    parser.add_argument("--ic-weight-scale", type=float, default=1.0, help="Initial physical weight multiplier for computational markers. Use 0 for source-only startup.")
    parser.add_argument("--pair-source", type=int, choices=[0, 1], default=1)
    parser.add_argument("--pair-source-weight-mode", type=int, choices=[0, 1], default=1, help="0 uses legacy BC_G weighting; 1 uses explicit pairSource_rate weighting.")
    parser.add_argument("--restart-path", type=Path, default=None, help="Previous PICOS++ output/HDF5 directory to load before normalization.")
    parser.add_argument("--restart-snapshot", type=int, default=-1, help="HDF5 snapshot index to load; -1 loads the latest numeric snapshot.")
    parser.add_argument("--restart-continue-time", action=argparse.BooleanOptionalAction, default=False, help="Continue the physical clock from the restart snapshot time.")
    parser.add_argument("--ic-type", type=int, default=1)
    parser.add_argument("--radius", type=float, default=0.05)
    parser.add_argument("--dtc", type=float, default=0.05)
    parser.add_argument("--mpis-for-fields", type=int, default=1)
    parser.add_argument("--mpi-ranks", type=int, default=2)
    parser.add_argument("--omp-threads", type=int, default=1)
    parser.add_argument("--fortran-omp-threads", type=int, default=None)
    parser.add_argument("--picos-omp-threads", type=int, default=None)
    parser.add_argument("--collisions", type=int, choices=[0, 1], default=1)
    parser.add_argument("--picos-run-mode", choices=["both", "nonrel", "rel"], default="both")
    parser.add_argument("--setup", action="store_true")
    parser.add_argument("--run-fortran", action="store_true")
    parser.add_argument("--run-picos", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = normalize_args(parser.parse_args())

    if args.setup:
        setup_cases(args)
    if args.run_fortran:
        if not fortran_case_path(args).is_file():
            setup_cases(args)
        run_fortran(args)
    if args.run_picos:
        if not (picos_input_dir(args) / f"input_file_{args.picos_tag_nonrel}.input").is_file():
            setup_cases(args)
        run_picos(args)
    plots: list[Path] = []
    if args.plot or not (args.setup or args.run_fortran or args.run_picos):
        plots = plot_comparison(args)
    summary = write_summary(args, plots)
    print(summary)
    for plot in plots:
        print(plot)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
