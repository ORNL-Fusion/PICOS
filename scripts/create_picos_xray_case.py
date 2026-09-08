#!/usr/bin/env python3
"""Create a PICOS input deck from an archived 2020 Proto-MPEX X-ray/ECH case.

This is intended to keep the benchmark setup reproducible. It converts the
archived namelist parameters and B(z) table into PICOS inputFiles entries:

  input_file_<tag>.input
  ions_properties_<tag>.ion
  xray_<case>_B_norm.txt
  xray_<case>_one.txt

The generated deck defaults to the PICOS 1D-2V guiding-center model
(`advanceParticleMethod=1`) and disables the field solve for direct comparison
to the archived X-ray cases. Poisson and 1D-3V Boris mode remain available
through explicit command-line flags.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Any

import numpy as np


E_CHARGE = 1.602176634e-19
EPS0 = 8.8541878128e-12
F_U = 1.66053906660e-27
F_C = 299792458.0


def parse_case_input(path: Path) -> dict[str, Any]:
    values: dict[str, Any] = {}
    pattern = re.compile(r"in%([A-Za-z0-9_]+)\s*=\s*([^!,]+)")
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("!", 1)[0].strip()
        if not line or line == "/" or "=" not in line:
            continue
        match = pattern.search(line)
        if not match:
            continue
        values[match.group(1)] = parse_value(match.group(2).strip().rstrip(","))
    return values


def parse_value(value: str) -> Any:
    value = value.strip().strip('"').strip("'")
    low = value.lower()
    if low in {".true.", "true"}:
        return True
    if low in {".false.", "false"}:
        return False
    try:
        if re.search(r"[.eEdD+-]", value):
            return float(value.replace("D", "E").replace("d", "e"))
        return int(value)
    except ValueError:
        return value


def find_bfield_file(case_dir: Path, xray_root: Path, metadata: dict[str, Any]) -> Path:
    name = str(metadata["BFieldFile"]).strip("/")
    candidates = [
        case_dir / name,
        case_dir / Path(name).name,
        xray_root / "BfieldData" / name,
        xray_root / "BfieldData" / Path(name).name,
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"Could not locate B-field file for {case_dir}")


def ion_skin_depth(ne_m3: float, mass_amu: float = 2.0, charge_state: float = 1.0) -> float:
    mass_kg = mass_amu * F_U
    omega_pi = math.sqrt(ne_m3 * (charge_state * E_CHARGE) ** 2 / (EPS0 * mass_kg))
    return F_C / omega_pi


def ion_gyroperiod(b_t: float, mass_amu: float = 2.0, charge_state: float = 1.0) -> float:
    mass_kg = mass_amu * F_U
    omega_ci = abs(charge_state) * E_CHARGE * b_t / mass_kg
    return 2.0 * math.pi / omega_ci


def enforce_mesh_multiple(nx: int, mpis_for_fields: int) -> int:
    multiple = math.lcm(2, mpis_for_fields)
    return max(multiple, nx - (nx % multiple))


def write_profile(path: Path, values: np.ndarray) -> None:
    np.savetxt(path, values, fmt="%.16e")


def build_input_text(
    metadata: dict[str, Any],
    tag: str,
    b_norm_file: str,
    one_file: str,
    args: argparse.Namespace,
    cv_b: float,
    dp: float,
    simulation_time_gyro: float,
    output_cadence_gyro: float,
    profile_points: int,
) -> str:
    sw_e = 1 if args.field_solve in {"poisson", "ohm"} else 0
    field_model = 1 if args.field_solve == "poisson" else 0
    rf_on = 1 if bool(metadata.get("iHeat", False)) else 0
    if args.force_rf_on:
        rf_on = 1
    if args.force_rf_off:
        rf_on = 0

    return f"""// PICOS input generated from 2020 X-ray benchmark {metadata.get('fileDescriptor', tag)}
// Source case keeps archived RF/B-field parameters but uses PICOS particle and field switches.
// =============================================================================
mpisForFields               {args.mpis_for_fields}
quietStart                  {args.quiet_start}
numberOfParticleSpecies     2
numberOfTracerSpecies       0
advanceParticleMethod       {args.advance_particle_method}

// Characteristic values:
// =============================================================================
CV_ne                       {float(metadata['ne0']):.16e}
CV_Te                       {float(metadata['Te0']):.16e}
CV_B                        {cv_b:.16e}
CV_Tpar                     {float(metadata['Ti0']):.16e}
CV_Tper                     {float(metadata['Ti0']):.16e}

// Simulation time:
// =============================================================================
DTc                         {args.dtc:.16e}
simulationTime              {simulation_time_gyro:.16e}

// Switches:
// =============================================================================
SW_EfieldSolve              {sw_e}
SW_fieldSolveModel          {field_model}
SW_BfieldSolve              0
SW_Collisions               {int(args.collisions)}
SW_RFheating                {rf_on}
SW_RFheatingIons            {args.rf_heat_ions}
SW_RFheatingElectrons       {args.rf_heat_electrons}
SW_relativisticElectrons    {args.relativistic_electrons}
SW_advancePos               1
SW_linearSolve              0

// Magnetic field initial conditions:
// =============================================================================
IC_uniformBfield            0
IC_BX                       {cv_b:.16e}
IC_BY                       0.0
IC_BZ                       0.0
IC_BX_NX                    {profile_points}
IC_BX_fileName              {b_norm_file}
IC_phiLeft                  {args.phi_left:.16e}
IC_phiRight                 {args.phi_right:.16e}
Poisson_BCModel             {args.poisson_bc_model}
Poisson_sheathCoefficient   {args.sheath_coefficient:.16e}

// Geometry:
// =============================================================================
dp                          {dp:.16e}
r1                          0.0
r2                          {args.radius:.16e}
LX_min                      {float(metadata['zmin']):.16e}
LX_max                      {float(metadata['zmax']):.16e}

// Electron initial conditions:
// =============================================================================
IC_ne                       {float(metadata['ne0']):.16e}
IC_Te                       {float(metadata['Te0']):.16e}
IC_Te_NX                    {profile_points}
IC_Te_fileName              {one_file}

// Ion RF operator:
// =============================================================================
RF_ion_Prf                      {args.rf_power:.16e}
RF_ion_n_harmonic               {int(metadata['n_harmonic'])}
RF_ion_freq                     {float(metadata['f_RF']):.16e}
RF_ion_x1                       {float(metadata['zRes1']):.16e}
RF_ion_x2                       {float(metadata['zRes2']):.16e}
RF_ion_t_ON                     0.0
RF_ion_t_OFF                    {simulation_time_gyro:.16e}
RF_ion_kpar                     {float(metadata['kpar']):.16e}
RF_ion_kper                     {float(metadata['kper']):.16e}
RF_ion_handedness               -1
RF_ion_EfieldMode               {args.rf_efield_mode}
RF_ion_EfieldAmplitude          {float(metadata['Ew']):.16e}
RF_ion_maxEnergyGainFraction    {args.rf_max_energy_gain_fraction:.16e}
RF_ion_maxParticleEnergy        {args.rf_max_particle_energy:.16e}
RF_ion_maxVelocityFractionC     {args.rf_max_velocity_fraction_c:.16e}
RF_ion_Prf_fileName             {one_file}
RF_ion_Prf_NS                   {profile_points}

// Electron RF/ECH operator:
// =============================================================================
RF_electron_Prf                      {args.rf_power:.16e}
RF_electron_n_harmonic               {int(metadata['n_harmonic'])}
RF_electron_freq                     {float(metadata['f_RF']):.16e}
RF_electron_x1                       {float(metadata['zRes1']):.16e}
RF_electron_x2                       {float(metadata['zRes2']):.16e}
RF_electron_t_ON                     0.0
RF_electron_t_OFF                    {simulation_time_gyro:.16e}
RF_electron_kpar                     {float(metadata['kpar']):.16e}
RF_electron_kper                     {float(metadata['kper']):.16e}
RF_electron_handedness               -1
RF_electron_EfieldMode               {args.rf_efield_mode}
RF_electron_EfieldAmplitude          {float(metadata['Ew']):.16e}
RF_electron_maxEnergyGainFraction    {args.rf_max_energy_gain_fraction:.16e}
RF_electron_maxParticleEnergy        {args.rf_max_particle_energy:.16e}
RF_electron_maxVelocityFractionC     {args.rf_max_velocity_fraction_c:.16e}
RF_electron_Prf_fileName             {one_file}
RF_electron_Prf_NS                   {profile_points}

// Output variables:
// =============================================================================
outputCadence               {output_cadence_gyro:.16e}
outputs_variables           {{X_p,V_p,a_p,BX_p,EX_p,BX_m,dBX_m,ddBX_m,n_m,Tpar_m,Tper_m,u_m,EX_m,Phi_m}}

// Data smoothing:
// =============================================================================
smoothingParameter          {args.smoothing:.16e}
filtersPerIterationFields   {args.field_filters}
filtersPerIterationIons     {args.ion_filters}
"""


def build_ions_text(metadata: dict[str, Any], one_file: str, profile_points: int, args: argparse.Namespace) -> str:
    ion_temp = float(metadata["Ti0"])
    electron_temp = float(metadata["Te0"])
    z_min = float(metadata["zmin"])
    z_max = float(metadata["zmax"])
    center = 0.5 * (z_min + z_max)
    sigma = 0.1 * (z_max - z_min)

    return f"""// PICOS species generated from 2020 X-ray benchmark {metadata.get('fileDescriptor', '')}
// =============================================================================
// Species 1: deuterium ions
// =============================================================================
SPECIES1                      1
NPC1                          {args.ion_particles_per_cell}
pctSupPartOutput1             {args.output_particle_percent:.16e}
Z1                            1
M1                            {float(metadata.get('Aion', 2.0)):.16e}

IC_type_1                     1
IC_Tper_1                     {ion_temp:.16e}
IC_Tper_fileName_1            {one_file}
IC_Tper_NX_1                  {profile_points}
IC_Tpar_1                     {ion_temp:.16e}
IC_Tpar_fileName_1            {one_file}
IC_Tpar_NX_1                  {profile_points}
IC_densityFraction_1          1.0
IC_densityFraction_fileName_1 {one_file}
IC_densityFraction_NX_1       {profile_points}

BC_type_1                     {args.boundary_type}
BC_T_1                        {ion_temp:.16e}
BC_E_1                        0.0
BC_eta_1                      0.7853981633974483
BC_mean_x_1                   {center:.16e}
BC_sigma_x_1                  {sigma:.16e}
BC_G_1                        0.0
BC_G_fileName_1               {one_file}
BC_G_NS_1                     {profile_points}

// =============================================================================
// Species 2: kinetic electrons
// =============================================================================
SPECIES2                      1
NPC2                          {args.electron_particles_per_cell}
pctSupPartOutput2             {args.output_particle_percent:.16e}
Z2                            -1
M2                            5.485799090441e-4

IC_type_2                     1
IC_Tper_2                     {electron_temp:.16e}
IC_Tper_fileName_2            {one_file}
IC_Tper_NX_2                  {profile_points}
IC_Tpar_2                     {electron_temp:.16e}
IC_Tpar_fileName_2            {one_file}
IC_Tpar_NX_2                  {profile_points}
IC_densityFraction_2          1.0
IC_densityFraction_fileName_2 {one_file}
IC_densityFraction_NX_2       {profile_points}

BC_type_2                     {args.boundary_type}
BC_T_2                        {electron_temp:.16e}
BC_E_2                        0.0
BC_eta_2                      0.7853981633974483
BC_mean_x_2                   {center:.16e}
BC_sigma_x_2                  {sigma:.16e}
BC_G_2                        0.0
BC_G_fileName_2               {one_file}
BC_G_NS_2                     {profile_points}
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xray-root", type=Path, default=Path("/Users/78k/Desktop/2020_07_13_XrayStudy"))
    parser.add_argument("--case", default="Case8")
    parser.add_argument("--tag", default=None)
    parser.add_argument("--input-dir", type=Path, default=Path("picosFILES/inputFiles"))
    parser.add_argument("--nx", type=int, default=200)
    parser.add_argument("--mpis-for-fields", type=int, default=2)
    parser.add_argument("--quiet-start", type=int, choices=[0, 1], default=1)
    parser.add_argument("--physical-time", type=float, default=None, help="Simulation time in seconds. Default uses archived Nsteps*dt.")
    parser.add_argument("--output-saves", type=int, default=5)
    parser.add_argument("--dtc", type=float, default=0.05)
    parser.add_argument("--field-solve", choices=["none", "ohm", "poisson"], default="none")
    parser.add_argument("--poisson-bc-model", type=int, default=1, help="0 Dirichlet, 1 periodic, 2 sheath")
    parser.add_argument("--sheath-coefficient", type=float, default=3.0)
    parser.add_argument("--phi-left", type=float, default=0.0)
    parser.add_argument("--phi-right", type=float, default=0.0)
    parser.add_argument("--advance-particle-method", type=int, default=1)
    parser.add_argument("--rf-efield-mode", type=int, default=1, help="0 absorbed-power balance, 1 fixed archived Ew")
    parser.add_argument("--rf-power", type=float, default=1.0)
    parser.add_argument("--rf-max-energy-gain-fraction", type=float, default=0.25)
    parser.add_argument("--rf-max-particle-energy", type=float, default=5000.0, help="Input in eV; 0 disables energy cap.")
    parser.add_argument("--rf-max-velocity-fraction-c", type=float, default=0.2, help="0 disables speed cap.")
    parser.add_argument("--collisions", type=int, default=1)
    parser.add_argument("--force-rf-on", action="store_true")
    parser.add_argument("--force-rf-off", action="store_true")
    parser.add_argument("--rf-heat-ions", type=int, choices=[0, 1], default=0)
    parser.add_argument("--rf-heat-electrons", type=int, choices=[0, 1], default=1)
    parser.add_argument("--relativistic-electrons", type=int, choices=[0, 1], default=1)
    parser.add_argument("--ion-particles-per-cell", type=int, default=4)
    parser.add_argument("--electron-particles-per-cell", type=int, default=16)
    parser.add_argument("--output-particle-percent", type=float, default=25.0)
    parser.add_argument("--boundary-type", type=int, default=3, help="3 periodic; use 4 for simple reinjection")
    parser.add_argument("--radius", type=float, default=0.05)
    parser.add_argument("--smoothing", type=float, default=1.0e-3)
    parser.add_argument("--field-filters", type=int, default=2)
    parser.add_argument("--ion-filters", type=int, default=2)
    args = parser.parse_args()

    xray_root = args.xray_root.expanduser().resolve()
    case_dir = xray_root / "OutputFiles" / "xp_Xray" / args.case
    metadata = parse_case_input(case_dir / "xp_Xray.in")
    bfield_path = find_bfield_file(case_dir, xray_root, metadata)

    tag = args.tag or f"xray_{args.case.lower()}"
    input_dir = args.input_dir.resolve()
    input_dir.mkdir(parents=True, exist_ok=True)

    nx = enforce_mesh_multiple(args.nx, args.mpis_for_fields)
    profile_points = nx + 2
    z_min = float(metadata["zmin"])
    z_max = float(metadata["zmax"])
    length = z_max - z_min

    b_data = np.loadtxt(bfield_path)
    z_raw = b_data[:, 0]
    b_raw = b_data[:, 1]
    order = np.argsort(z_raw)
    z_raw = z_raw[order]
    b_raw = b_raw[order]

    dx_profile = length / nx
    z_profile = z_min - 0.5 * dx_profile + dx_profile * np.arange(profile_points)
    b_interp = np.interp(z_profile, z_raw, b_raw)
    cv_b = float(np.nanmax(np.abs(b_interp)))
    b_norm = b_interp / cv_b
    one = np.ones(profile_points)

    physical_time = args.physical_time
    if physical_time is None:
        physical_time = float(metadata["Nsteps"]) * float(metadata["dt"])
    tgyro = ion_gyroperiod(cv_b, mass_amu=float(metadata.get("Aion", 2.0)), charge_state=float(metadata.get("Zion", 1.0)))
    simulation_time_gyro = physical_time / tgyro
    output_cadence_gyro = simulation_time_gyro / max(args.output_saves, 1)
    dp = (length / nx) / ion_skin_depth(float(metadata["ne0"]), mass_amu=float(metadata.get("Aion", 2.0)), charge_state=float(metadata.get("Zion", 1.0)))

    b_norm_file = f"{tag}_B_norm.txt"
    one_file = f"{tag}_one.txt"
    write_profile(input_dir / b_norm_file, b_norm)
    write_profile(input_dir / one_file, one)

    input_text = build_input_text(
        metadata=metadata,
        tag=tag,
        b_norm_file=b_norm_file,
        one_file=one_file,
        args=args,
        cv_b=cv_b,
        dp=dp,
        simulation_time_gyro=simulation_time_gyro,
        output_cadence_gyro=output_cadence_gyro,
        profile_points=profile_points,
    )
    ions_text = build_ions_text(metadata, one_file, profile_points, args)

    (input_dir / f"input_file_{tag}.input").write_text(input_text)
    (input_dir / f"ions_properties_{tag}.ion").write_text(ions_text)

    print(f"Wrote {input_dir / f'input_file_{tag}.input'}")
    print(f"Wrote {input_dir / f'ions_properties_{tag}.ion'}")
    print(f"Wrote {input_dir / b_norm_file}")
    print(f"Wrote {input_dir / one_file}")
    print(f"PICOS tag: {tag}")
    print(f"Approximate NX requested/generated: {args.nx}/{nx}")
    print(f"Archived physical time represented: {physical_time:.6e} s = {simulation_time_gyro:.6e} ion gyroperiods")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
