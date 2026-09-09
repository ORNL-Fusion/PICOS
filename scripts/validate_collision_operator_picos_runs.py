#!/usr/bin/env python3
"""End-to-end PICOS++ collision validations.

This driver creates small collision-only PICOS++ input decks, runs the real
`xpicos` executable, reads particle HDF5 outputs with `h5dump`, and checks:

* thermal stability of a single kinetic species,
* hot/cold e-e equilibration,
* hot-electron/cold-ion e-i equilibration,
* hot/cold i-i equilibration,
* slowing down of a hot trace ion population,
* particle/charge conservation in all closed periodic cases.

The cases intentionally disable field solve, RF, sources, and axial motion so
the only active physics in the normal timestep loop is the collision operator.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


E_CHARGE = 1.602176634e-19
M_E = 9.1093837015e-31
M_U = 1.66053906660e-27
EPS0 = 8.8541878128e-12
C_LIGHT = 299792458.0


@dataclass(frozen=True)
class SpeciesDeck:
    label: str
    z: float
    mass_amu: float
    density_fraction: float
    tpar_ev: float
    tper_ev: float
    particles_per_cell: int


@dataclass(frozen=True)
class CaseDeck:
    tag: str
    description: str
    species: tuple[SpeciesDeck, ...]
    density_m3: float
    electron_fluid_te_ev: float
    b_t: float
    lx_min: float
    lx_max: float
    radius_m: float
    nx: int
    profile_n: int
    dtc: float
    physical_time_s: float
    output_intervals: int


def default_cases(particles_per_cell: int) -> list[CaseDeck]:
    density = 1.0e20
    b_t = 0.01
    nx = 16
    profile_n = 34
    hot_ppc = max(8, particles_per_cell)
    cold_ppc = max(8, particles_per_cell)
    trace_ppc = max(8, particles_per_cell // 2)
    return [
        CaseDeck(
            tag="coll_thermal_stability_ee",
            description="single kinetic electron species, Maxwellian self-collision stability",
            species=(
                SpeciesDeck("e thermal", -1.0, M_E / M_U, 1.0, 60.0, 60.0, hot_ppc),
            ),
            density_m3=density,
            electron_fluid_te_ev=60.0,
            b_t=b_t,
            lx_min=-0.01,
            lx_max=0.01,
            radius_m=0.03,
            nx=nx,
            profile_n=profile_n,
            dtc=800.0,
            physical_time_s=1.0e-7,
            output_intervals=10,
        ),
        CaseDeck(
            tag="coll_equilibration_ee",
            description="two kinetic electron populations equilibrating by e-e collisions",
            species=(
                SpeciesDeck("e cold", -1.0, M_E / M_U, 0.5, 20.0, 20.0, cold_ppc),
                SpeciesDeck("e hot", -1.0, M_E / M_U, 0.5, 180.0, 180.0, hot_ppc),
            ),
            density_m3=density,
            electron_fluid_te_ev=100.0,
            b_t=b_t,
            lx_min=-0.01,
            lx_max=0.01,
            radius_m=0.03,
            nx=nx,
            profile_n=profile_n,
            dtc=800.0,
            physical_time_s=1.5e-7,
            output_intervals=15,
        ),
        CaseDeck(
            tag="coll_equilibration_ei",
            description="accelerated high-density cold-D+/hot-electron e-i equilibration check",
            species=(
                SpeciesDeck("D+ cold", 1.0, 2.01410177812, 1.0, 5.0, 5.0, cold_ppc),
                SpeciesDeck("e hot", -1.0, M_E / M_U, 1.0, 180.0, 180.0, hot_ppc),
            ),
            density_m3=6.0e22,
            electron_fluid_te_ev=180.0,
            b_t=b_t,
            lx_min=-0.01,
            lx_max=0.01,
            radius_m=0.03,
            nx=nx,
            profile_n=profile_n,
            dtc=800.0,
            physical_time_s=2.0e-7,
            output_intervals=16,
        ),
        CaseDeck(
            tag="coll_equilibration_ii",
            description="hot/cold D+ ion populations equilibrating by i-i collisions",
            species=(
                SpeciesDeck("D+ cold", 1.0, 2.01410177812, 0.5, 10.0, 10.0, cold_ppc),
                SpeciesDeck("D+ hot", 1.0, 2.01410177812, 0.5, 180.0, 180.0, hot_ppc),
            ),
            density_m3=density,
            electron_fluid_te_ev=95.0,
            b_t=b_t,
            lx_min=-0.01,
            lx_max=0.01,
            radius_m=0.03,
            nx=nx,
            profile_n=profile_n,
            dtc=800.0,
            physical_time_s=6.0e-6,
            output_intervals=20,
        ),
        CaseDeck(
            tag="coll_slowing_down_trace_ion",
            description="hot trace D+ slowing down on cold D+ background",
            species=(
                SpeciesDeck("D+ background", 1.0, 2.01410177812, 1.0, 8.0, 8.0, cold_ppc),
                SpeciesDeck("D+ trace hot", 1.0, 2.01410177812, 0.05, 500.0, 500.0, trace_ppc),
            ),
            density_m3=density,
            electron_fluid_te_ev=8.0,
            b_t=b_t,
            lx_min=-0.01,
            lx_max=0.01,
            radius_m=0.03,
            nx=nx,
            profile_n=profile_n,
            dtc=800.0,
            physical_time_s=5.0e-6,
            output_intervals=20,
        ),
    ]


def reference_species(case: CaseDeck) -> SpeciesDeck:
    for species in case.species:
        if species.z > 0.0:
            return species
    return case.species[0]


def gyroperiod_s(species: SpeciesDeck, b_t: float) -> float:
    mass = species.mass_amu * M_U
    return 2.0 * math.pi * mass / (abs(species.z) * E_CHARGE * b_t)


def skin_depth_m(species: SpeciesDeck, density_m3: float) -> float:
    mass = species.mass_amu * M_U
    n = max(density_m3 * species.density_fraction, 1.0)
    wp = math.sqrt(n * (species.z * E_CHARGE) ** 2 / (EPS0 * mass))
    return C_LIGHT / wp


def write_flat_profile(path: Path, n_values: int, value: float = 1.0) -> None:
    path.write_text("\n".join(f"{value:.16e}" for _ in range(n_values)) + "\n")


def input_text(case: CaseDeck) -> str:
    ref = reference_species(case)
    ref_gyro = gyroperiod_s(ref, case.b_t)
    ref_skin = skin_depth_m(ref, case.density_m3)
    sim_gyroperiods = case.physical_time_s / ref_gyro
    output_gyroperiods = sim_gyroperiods / max(case.output_intervals, 1)
    dx = (case.lx_max - case.lx_min) / case.nx
    dp = dx / ref_skin
    profile = f"{case.tag}_one.txt"
    return f"""// PICOS++ collision validation deck: {case.description}
// Generated by scripts/validate_collision_operator_picos_runs.py.
// Physics switches intentionally isolate Coulomb collisions:
//   no field solve, no RF, no pair source, no axial motion, periodic particles.
// Units: temperatures in eV, density in m^-3, magnetic field in T.
// =============================================================================
mpisForFields               1
quietStart                  1
IC_velocityDistributionModel 0
IC_randomSeed              271828
numberOfParticleSpecies     {len(case.species)}
numberOfTracerSpecies       0
advanceParticleMethod       1

// Characteristic values:
// =============================================================================
CV_ne                       {case.density_m3:.16e}
CV_Te                       {case.electron_fluid_te_ev:.16e}
CV_B                        {case.b_t:.16e}
CV_Tpar                     {case.species[0].tpar_ev:.16e}
CV_Tper                     {case.species[0].tper_ev:.16e}

// Simulation time:
// =============================================================================
DTc                         {case.dtc:.16e}
simulationTime              {sim_gyroperiods:.16e}

// Switches:
// =============================================================================
SW_EfieldSolve              0
SW_fieldSolveModel          0
SW_BfieldSolve              0
SW_Collisions               1
CollOperType                2
SW_collisionConservationProjection 1
collisionRandomSeed        314159
SW_RFheating                0
SW_RFheatingIons            0
SW_RFheatingElectrons       0
SW_pairSource               0
SW_relativisticElectrons    0
SW_advancePos               0
SW_linearSolve              0

// Magnetic field initial conditions:
// =============================================================================
IC_uniformBfield            1
IC_BX                       {case.b_t:.16e}
IC_BY                       0.0
IC_BZ                       0.0
IC_BX_NX                    {case.profile_n}
IC_BX_fileName              {profile}
IC_phiLeft                  0.0
IC_phiRight                 0.0
Poisson_BCModel             1
Poisson_sheathCoefficient   3.0
ReformulatedPoisson_lambda  -1.0
ReformulatedPoisson_quasiNeutral 0

// Geometry:
// =============================================================================
dp                          {dp:.16e}
r1                          0.0
r2                          {case.radius_m:.16e}
LX_min                      {case.lx_min:.16e}
LX_max                      {case.lx_max:.16e}

// Electron fluid/profile initial conditions:
// =============================================================================
IC_ne                       {case.density_m3:.16e}
IC_Te                       {case.electron_fluid_te_ev:.16e}
IC_Te_NX                    {case.profile_n}
IC_Te_fileName              {profile}

// Coupled electron-ion source, disabled for closed collision validation:
// =============================================================================
pairSource_ionSpecies       1
pairSource_electronSpecies  2
pairSource_rate             0.0
pairSource_mean_x           0.0
pairSource_sigma_x          0.001
pairSource_Ti_birth         {case.electron_fluid_te_ev:.16e}
pairSource_Te_birth         {case.electron_fluid_te_ev:.16e}
pairSource_Ei_birth         0.0
pairSource_Ee_birth         0.0
pairSource_eta_i            0.0
pairSource_eta_e            0.0
pairSource_positionMode     0
pairSource_fileName         {profile}
pairSource_NS               {case.profile_n}
pairSource_maxParticleWeight 1000

// RF operator, disabled for closed collision validation:
// =============================================================================
RF_ion_Prf                      0.0
RF_ion_n_harmonic               1
RF_ion_freq                     1.0
RF_ion_x1                       {case.lx_min:.16e}
RF_ion_x2                       {case.lx_max:.16e}
RF_ion_t_ON                     0.0
RF_ion_t_OFF                    {sim_gyroperiods:.16e}
RF_ion_kpar                     0.0
RF_ion_kper                     0.0
RF_ion_handedness               -1
RF_ion_EfieldMode               1
RF_ion_EfieldAmplitude          0.0
RF_ion_maxEnergyGainFraction    0.0
RF_ion_maxParticleEnergy        0.0
RF_ion_maxVelocityFractionC     0.0
RF_ion_Prf_fileName             {profile}
RF_ion_Prf_NS                   {case.profile_n}

RF_electron_Prf                      0.0
RF_electron_n_harmonic               1
RF_electron_freq                     1.0
RF_electron_x1                       {case.lx_min:.16e}
RF_electron_x2                       {case.lx_max:.16e}
RF_electron_t_ON                     0.0
RF_electron_t_OFF                    {sim_gyroperiods:.16e}
RF_electron_kpar                     0.0
RF_electron_kper                     0.0
RF_electron_handedness               -1
RF_electron_EfieldMode               1
RF_electron_EfieldAmplitude          0.0
RF_electron_maxEnergyGainFraction    0.0
RF_electron_maxParticleEnergy        0.0
RF_electron_maxVelocityFractionC     0.0
RF_electron_Prf_fileName             {profile}
RF_electron_Prf_NS                   {case.profile_n}

// Output variables:
// =============================================================================
outputCadence               {output_gyroperiods:.16e}
outputs_variables           {{X_p,V_p,a_p,BX_p,n_m,Tpar_m,Tper_m,u_m}}

// Data smoothing:
// =============================================================================
smoothingParameter          0.0
filtersPerIterationFields   0
filtersPerIterationIons     0
"""


def ions_text(case: CaseDeck) -> str:
    profile = f"{case.tag}_one.txt"
    blocks = [
        f"// PICOS++ collision validation species: {case.description}",
        "// Generated by scripts/validate_collision_operator_picos_runs.py.",
        "// =============================================================================",
    ]
    for idx, species in enumerate(case.species, start=1):
        blocks.append(
            f"""
// Species {idx}: {species.label}
// =============================================================================
SPECIES{idx}                      1
NPC{idx}                          {species.particles_per_cell}
pctSupPartOutput{idx}             100
Z{idx}                            {species.z:.16e}
M{idx}                            {species.mass_amu:.16e}

IC_type_{idx}                     1
IC_Tper_{idx}                     {species.tper_ev:.16e}
IC_Tper_fileName_{idx}            {profile}
IC_Tper_NX_{idx}                  {case.profile_n}
IC_Tpar_{idx}                     {species.tpar_ev:.16e}
IC_Tpar_fileName_{idx}            {profile}
IC_Tpar_NX_{idx}                  {case.profile_n}
IC_densityFraction_{idx}          {species.density_fraction:.16e}
IC_densityFraction_fileName_{idx} {profile}
IC_densityFraction_NX_{idx}       {case.profile_n}

BC_type_{idx}                     3
BC_T_{idx}                        {species.tpar_ev:.16e}
BC_E_{idx}                        0.0
BC_eta_{idx}                      0.0
BC_mean_x_{idx}                   0.0
BC_sigma_x_{idx}                  0.001
BC_G_{idx}                        0.0
BC_G_fileName_{idx}               {profile}
BC_G_NS_{idx}                     {case.profile_n}
""".strip()
        )
    return "\n\n".join(blocks) + "\n"


def h5dump_text(path: Path, dataset: str | None = None, header_only: bool = False, names: bool = False) -> str:
    cmd = ["h5dump"]
    if names:
        cmd.append("-n")
    if header_only:
        cmd.append("-H")
    if dataset is not None:
        cmd.extend(["-d", dataset])
    cmd.append(str(path))
    return subprocess.check_output(cmd, text=True)


def hdf5_steps(path: Path) -> list[str]:
    text = h5dump_text(path, names=True)
    steps = sorted({match.group(1) for match in re.finditer(r"group\s+/([0-9]+)(?:\s|$)", text)}, key=int)
    if not steps:
        raise RuntimeError(f"No output-step groups found in {path}")
    return steps


def read_hdf5_dataset(path: Path, dataset: str) -> np.ndarray:
    header = h5dump_text(path, dataset=dataset, header_only=True)
    shape_match = re.search(r"DATASPACE\s+SIMPLE\s+\{\s+\(\s*([0-9,\s]+)\)", header)
    if not shape_match:
        raise RuntimeError(f"Could not parse HDF5 shape for {dataset} in {path}")
    shape = tuple(int(item.strip()) for item in shape_match.group(1).split(",") if item.strip())
    if "H5T_IEEE_F64LE" in header:
        dtype = "<f8"
    elif "H5T_IEEE_F32LE" in header:
        dtype = "<f4"
    elif "H5T_STD_I32LE" in header:
        dtype = "<i4"
    else:
        raise RuntimeError(f"Unsupported HDF5 datatype for {dataset} in {path}")

    with tempfile.NamedTemporaryFile(prefix="picos_collision_h5_", suffix=".bin") as payload, tempfile.NamedTemporaryFile(
        prefix="picos_collision_h5_", suffix=".ddl"
    ) as ddl:
        subprocess.check_call(["h5dump", "-d", dataset, "-o", payload.name, "-b", "LE", "-O", ddl.name, str(path)])
        data = np.fromfile(payload.name, dtype=dtype)
    return data.reshape(shape)


def read_scalar(path: Path, dataset: str) -> float:
    data = read_hdf5_dataset(path, dataset)
    return float(np.ravel(data)[0])


def particle_files(hdf5_dir: Path) -> list[Path]:
    files = sorted(hdf5_dir.glob("PARTICLES_FILE_*.h5"), key=lambda item: int(item.stem.rsplit("_", 1)[1]))
    if not files:
        raise RuntimeError(f"No PARTICLES_FILE_*.h5 files found in {hdf5_dir}")
    return files


def velocity_components(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2:
        raise RuntimeError(f"Expected 2D V_p dataset, got {values.shape}")
    if values.shape[0] in (2, 3):
        return values
    if values.shape[1] in (2, 3):
        return values.T
    raise RuntimeError(f"Cannot infer velocity axis for V_p shape {values.shape}")


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(mask):
        return math.nan
    return float(np.average(values[mask], weights=weights[mask]))


def species_stats(hdf5_dir: Path, step: str, species_index: int) -> dict[str, float]:
    main = hdf5_dir / "main.h5"
    species_path = f"/ions/spp_{species_index}"
    mass = read_scalar(main, f"{species_path}/M")
    z = read_scalar(main, f"{species_path}/Z")
    ncp = read_scalar(main, f"{species_path}/NCP")

    velocities = []
    weights = []
    finite = True
    for path in particle_files(hdf5_dir):
        dataset = f"/{step}/ions/spp_{species_index}"
        try:
            velocities.append(velocity_components(read_hdf5_dataset(path, f"{dataset}/V_p")))
            weights.append(np.ravel(read_hdf5_dataset(path, f"{dataset}/a_p")).astype(float))
        except subprocess.CalledProcessError:
            continue
    if not velocities:
        raise RuntimeError(f"No particle velocity output for spp_{species_index} step {step} in {hdf5_dir}")

    v = np.concatenate(velocities, axis=1)
    a = np.concatenate(weights)
    finite = bool(np.isfinite(v).all() and np.isfinite(a).all())
    real_weight = a * ncp
    vpar = v[0]
    if v.shape[0] == 2:
        vperp2 = v[1] * v[1]
        speed2 = vpar * vpar + vperp2
    else:
        vperp2 = v[1] * v[1] + v[2] * v[2]
        speed2 = vpar * vpar + vperp2

    mean_vpar = weighted_mean(vpar, real_weight)
    tpar = mass * weighted_mean((vpar - mean_vpar) ** 2, real_weight) / E_CHARGE
    tper = 0.5 * mass * weighted_mean(vperp2, real_weight) / E_CHARGE
    scalar_t = (tpar + 2.0 * tper) / 3.0
    mean_energy_ev = 0.5 * mass * weighted_mean(speed2, real_weight) / E_CHARGE
    total_real = float(np.sum(real_weight))
    total_charge_number = float(z * total_real)
    total_energy_j = float(np.sum(0.5 * mass * speed2 * real_weight))
    return {
        "mass_kg": mass,
        "z": z,
        "ncp": ncp,
        "finite": finite,
        "super_particles": float(a.size),
        "real_particles": total_real,
        "charge_number": total_charge_number,
        "total_energy_j": total_energy_j,
        "mean_vpar_m_s": mean_vpar,
        "tpar_ev": tpar,
        "tper_ev": tper,
        "scalar_t_ev": scalar_t,
        "mean_energy_ev": mean_energy_ev,
    }


def read_case_trace(repo_root: Path, tag: str, n_species: int) -> list[dict[str, Any]]:
    hdf5_dir = repo_root / "picosFILES" / "outputFiles" / tag / "HDF5"
    first_particle = particle_files(hdf5_dir)[0]
    steps = hdf5_steps(first_particle)
    trace: list[dict[str, Any]] = []
    for step in steps:
        time_s = read_scalar(first_particle, f"/{step}/time")
        row: dict[str, Any] = {"step": int(step), "time_s": time_s, "species": []}
        for idx in range(1, n_species + 1):
            row["species"].append(species_stats(hdf5_dir, step, idx))
        trace.append(row)
    return trace


def write_trace_csv(path: Path, trace: list[dict[str, Any]]) -> None:
    max_species = max(len(row["species"]) for row in trace)
    fields = ["step", "time_s"]
    for idx in range(1, max_species + 1):
        for key in (
            "scalar_t_ev",
            "tpar_ev",
            "tper_ev",
            "mean_energy_ev",
            "total_energy_j",
            "real_particles",
            "charge_number",
            "mean_vpar_m_s",
        ):
            fields.append(f"spp_{idx}_{key}")

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in trace:
            out = {"step": row["step"], "time_s": f"{row['time_s']:.16e}"}
            for idx, species in enumerate(row["species"], start=1):
                for key in fields:
                    prefix = f"spp_{idx}_"
                    if key.startswith(prefix):
                        raw_key = key[len(prefix) :]
                        out[key] = f"{species[raw_key]:.16e}" if isinstance(species[raw_key], float) else species[raw_key]
            writer.writerow(out)


def total_for(trace_row: dict[str, Any], key: str) -> float:
    return float(sum(species[key] for species in trace_row["species"]))


def conservation_metrics(trace: list[dict[str, Any]]) -> dict[str, float | bool]:
    first = trace[0]
    rel_particle_drifts = []
    for idx, species0 in enumerate(first["species"]):
        initial = abs(species0["real_particles"])
        values = [row["species"][idx]["real_particles"] for row in trace]
        denom = max(initial, 1.0)
        rel_particle_drifts.append(max(abs(value - values[0]) / denom for value in values))

    charge0 = total_for(first, "charge_number")
    abs_charge0 = sum(abs(species["charge_number"]) for species in first["species"])
    charge_abs_drift = max(abs(total_for(row, "charge_number") - charge0) for row in trace)
    charge_rel_drift = charge_abs_drift / max(abs_charge0, 1.0)

    energy0 = total_for(first, "total_energy_j")
    energy_rel_drift = max(abs(total_for(row, "total_energy_j") - energy0) for row in trace) / max(abs(energy0), 1.0e-300)
    all_finite = all(bool(species["finite"]) for row in trace for species in row["species"])
    return {
        "all_finite": all_finite,
        "max_species_particle_rel_drift": float(max(rel_particle_drifts)),
        "total_charge_rel_drift": float(charge_rel_drift),
        "total_energy_rel_drift": float(energy_rel_drift),
    }


def gap_reduction(trace: list[dict[str, Any]], idx_a: int = 0, idx_b: int = 1) -> dict[str, float]:
    initial_a = trace[0]["species"][idx_a]["scalar_t_ev"]
    initial_b = trace[0]["species"][idx_b]["scalar_t_ev"]
    final_a = trace[-1]["species"][idx_a]["scalar_t_ev"]
    final_b = trace[-1]["species"][idx_b]["scalar_t_ev"]
    initial_gap = abs(initial_b - initial_a)
    final_gap = abs(final_b - final_a)
    return {
        "initial_t_a_ev": float(initial_a),
        "initial_t_b_ev": float(initial_b),
        "final_t_a_ev": float(final_a),
        "final_t_b_ev": float(final_b),
        "initial_gap_ev": float(initial_gap),
        "final_gap_ev": float(final_gap),
        "gap_ratio": float(final_gap / max(initial_gap, 1.0e-300)),
    }


def fit_slowing_down(trace: list[dict[str, Any]], hot_index: int, background_temperature_ev: float) -> dict[str, float]:
    times = np.asarray([row["time_s"] for row in trace], dtype=float)
    energies = np.asarray([row["species"][hot_index]["mean_energy_ev"] for row in trace], dtype=float)
    excess = np.maximum(energies - 1.5 * background_temperature_ev, 1.0e-12)
    mask = np.isfinite(times) & np.isfinite(excess) & (excess > 0.0)
    if np.count_nonzero(mask) < 3:
        return {"tau_fit_s": math.nan, "decay_rate_s_1": math.nan, "energy_drop_fraction": math.nan}
    slope, intercept = np.polyfit(times[mask], np.log(excess[mask]), 1)
    tau = -1.0 / slope if slope < 0.0 else math.inf
    drop = (excess[0] - excess[-1]) / max(excess[0], 1.0e-300)
    return {
        "tau_fit_s": float(tau),
        "decay_rate_s_1": float(-slope),
        "log_excess_intercept": float(intercept),
        "energy_drop_fraction": float(drop),
    }


def analyze_case(case: CaseDeck, trace: list[dict[str, Any]]) -> dict[str, Any]:
    cons = conservation_metrics(trace)
    summary: dict[str, Any] = {
        "description": case.description,
        "conservation": cons,
        "initial": [species["scalar_t_ev"] for species in trace[0]["species"]],
        "final": [species["scalar_t_ev"] for species in trace[-1]["species"]],
    }

    common_pass = bool(
        cons["all_finite"]
        and cons["max_species_particle_rel_drift"] < 1.0e-10
        and cons["total_charge_rel_drift"] < 1.0e-10
    )

    if case.tag == "coll_thermal_stability_ee":
        initial_t = trace[0]["species"][0]["scalar_t_ev"]
        final_t = trace[-1]["species"][0]["scalar_t_ev"]
        anisotropy = abs(trace[-1]["species"][0]["tpar_ev"] - trace[-1]["species"][0]["tper_ev"]) / max(final_t, 1.0e-300)
        rel_t_drift = abs(final_t - initial_t) / max(initial_t, 1.0e-300)
        summary["thermal_stability"] = {
            "relative_temperature_drift": float(rel_t_drift),
            "final_parallel_perp_anisotropy": float(anisotropy),
        }
        summary["passed"] = common_pass and rel_t_drift < 0.30 and anisotropy < 0.45
    elif case.tag in {"coll_equilibration_ee", "coll_equilibration_ei", "coll_equilibration_ii"}:
        gap = gap_reduction(trace)
        cold_warmed = gap["final_t_a_ev"] > 0.90 * gap["initial_t_a_ev"]
        hot_cooled = gap["final_t_b_ev"] < 1.10 * gap["initial_t_b_ev"]
        summary["equilibration"] = gap
        summary["passed"] = common_pass and gap["gap_ratio"] < 0.85 and cold_warmed and hot_cooled
    elif case.tag == "coll_slowing_down_trace_ion":
        fit = fit_slowing_down(trace, hot_index=1, background_temperature_ev=8.0)
        summary["slowing_down"] = fit
        summary["passed"] = common_pass and math.isfinite(fit["tau_fit_s"]) and fit["energy_drop_fraction"] > 0.10
    else:
        summary["passed"] = common_pass
    return summary


def stage_case(repo_root: Path, case: CaseDeck) -> None:
    input_dir = repo_root / "picosFILES" / "inputFiles"
    input_dir.mkdir(parents=True, exist_ok=True)
    write_flat_profile(input_dir / f"{case.tag}_one.txt", case.profile_n, 1.0)
    (input_dir / f"input_file_{case.tag}.input").write_text(input_text(case))
    (input_dir / f"ions_properties_{case.tag}.ion").write_text(ions_text(case))


def run_case(repo_root: Path, binary: Path, case: CaseDeck, mpi_ranks: int, keep_existing: bool, out_dir: Path) -> None:
    output_root = repo_root / "picosFILES" / "outputFiles"
    output_path = output_root / case.tag
    if output_path.exists() and not keep_existing:
        shutil.rmtree(output_path)

    cmd = ["mpirun", "-np", str(mpi_ranks), str(binary), "1-D", "outputFiles", case.tag]
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    result = subprocess.run(cmd, cwd=repo_root / "picosFILES", env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log_path = out_dir / f"{case.tag}.log"
    log_path.write_text(result.stdout)
    if result.returncode != 0:
        tail = "\n".join(result.stdout.splitlines()[-40:])
        raise RuntimeError(f"PICOS++ run failed for {case.tag} with exit code {result.returncode}.\n{tail}")


def make_plots(out_dir: Path, traces: dict[str, list[dict[str, Any]]]) -> list[Path]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from paper_plot_style import apply_paper_figure_style

        apply_paper_figure_style()
    except Exception:
        return []

    paths: list[Path] = []
    title_labels = {
        "coll_thermal_stability_ee": "e-e Stability",
        "coll_equilibration_ee": "e-e Relaxation",
        "coll_equilibration_ei": "e-i Relaxation",
        "coll_equilibration_ii": "i-i Relaxation",
        "coll_slowing_down_trace_ion": "Trace-Ion Slowing",
    }
    for tag, trace in traces.items():
        times = np.asarray([row["time_s"] for row in trace], dtype=float)
        fig, axes = plt.subplots(1, 2, figsize=(18, 6.5), constrained_layout=True)
        for idx in range(len(trace[0]["species"])):
            scalar_t = [row["species"][idx]["scalar_t_ev"] for row in trace]
            tpar = [row["species"][idx]["tpar_ev"] for row in trace]
            tper = [row["species"][idx]["tper_ev"] for row in trace]
            axes[0].plot(times, scalar_t, marker="o", ms=3, label=f"spp_{idx + 1} T")
            axes[0].plot(times, tpar, ls="--", alpha=0.45, label=f"spp_{idx + 1} Tpar")
            axes[0].plot(times, tper, ls=":", alpha=0.55, label=f"spp_{idx + 1} Tper")
        axes[0].set_xlabel("time [s]")
        axes[0].set_ylabel("temperature [eV]")
        axes[0].set_title(title_labels.get(tag, tag))
        axes[0].grid(True, alpha=0.25)
        axes[0].legend(fontsize=24)

        particles0 = np.asarray([sum(species["real_particles"] for species in row["species"]) for row in trace], dtype=float)
        charge0 = np.asarray([sum(species["charge_number"] for species in row["species"]) for row in trace], dtype=float)
        energy0 = np.asarray([sum(species["total_energy_j"] for species in row["species"]) for row in trace], dtype=float)
        axes[1].plot(times, (particles0 - particles0[0]) / max(abs(particles0[0]), 1.0), label="real particles")
        axes[1].plot(times, (charge0 - charge0[0]) / max(np.max(np.abs(charge0)), 1.0), label="net charge")
        axes[1].plot(times, (energy0 - energy0[0]) / max(abs(energy0[0]), 1.0e-300), label="kinetic energy")
        axes[1].set_xlabel("time [s]")
        axes[1].set_ylabel("relative change")
        axes[1].set_title("closed-case conservation")
        axes[1].grid(True, alpha=0.25)
        axes[1].legend(fontsize=24)
        path = out_dir / f"{tag}.png"
        fig.savefig(path, dpi=180, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)

    fig, axes = plt.subplots(2, 2, figsize=(18, 12), constrained_layout=True)
    selected = ["coll_equilibration_ee", "coll_equilibration_ei", "coll_equilibration_ii", "coll_slowing_down_trace_ion"]
    for ax, tag in zip(axes.ravel(), selected):
        if tag not in traces:
            ax.axis("off")
            continue
        trace = traces[tag]
        times = np.asarray([row["time_s"] for row in trace], dtype=float)
        if tag == "coll_slowing_down_trace_ion":
            values = [row["species"][1]["mean_energy_ev"] for row in trace]
            ax.plot(times, values, marker="o", ms=3)
            ax.set_ylabel("trace mean energy [eV]")
        else:
            for idx in range(len(trace[0]["species"])):
                values = [row["species"][idx]["scalar_t_ev"] for row in trace]
                ax.plot(times, values, marker="o", ms=3, label=f"spp_{idx + 1}")
            ax.set_ylabel("temperature [eV]")
            ax.legend(fontsize=24)
        ax.set_title(title_labels.get(tag, tag))
        ax.set_xlabel("time [s]")
        ax.grid(True, alpha=0.25)
    path = out_dir / "picos_collision_run_validation.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--picos-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--binary", type=Path, default=None, help="Path to xpicos. Default: build/picosFILES/src/xpicos.")
    parser.add_argument("--output-dir", type=Path, default=Path("validation/collision_operator_picos_runs"))
    parser.add_argument("--mpi-ranks", type=int, default=2)
    parser.add_argument("--particles-per-cell", type=int, default=64)
    parser.add_argument("--skip-run", action="store_true", help="Analyze existing outputFiles/<tag>/HDF5 directories.")
    parser.add_argument("--keep-existing", action="store_true", help="Do not remove existing generated output directories before running.")
    parser.add_argument("--case", action="append", dest="cases", help="Run/analyze only this case tag. Can be repeated.")
    args = parser.parse_args()

    repo_root = args.picos_root.resolve()
    binary = (repo_root / "build" / "picosFILES" / "src" / "xpicos") if args.binary is None else args.binary.resolve()
    out_dir = args.output_dir if args.output_dir.is_absolute() else repo_root / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    cases = default_cases(args.particles_per_cell)
    if args.cases:
        requested = set(args.cases)
        cases = [case for case in cases if case.tag in requested]
        missing = requested - {case.tag for case in cases}
        if missing:
            raise SystemExit(f"Unknown case tag(s): {', '.join(sorted(missing))}")
    if not cases:
        raise SystemExit("No cases selected")

    if args.mpi_ranks % 2:
        raise SystemExit("--mpi-ranks must be even for PICOS++")
    if not args.skip_run and not binary.is_file():
        raise SystemExit(f"xpicos binary not found: {binary}")

    traces: dict[str, list[dict[str, Any]]] = {}
    summaries: dict[str, Any] = {}
    for case in cases:
        stage_case(repo_root, case)
        if not args.skip_run:
            print(f"Running {case.tag} ...", flush=True)
            run_case(repo_root, binary, case, args.mpi_ranks, args.keep_existing, out_dir)
        trace = read_case_trace(repo_root, case.tag, len(case.species))
        traces[case.tag] = trace
        write_trace_csv(out_dir / f"{case.tag}.csv", trace)
        summaries[case.tag] = analyze_case(case, trace)

    plot_paths = make_plots(out_dir, traces)
    overall_passed = all(bool(summary.get("passed", False)) for summary in summaries.values())
    summary = {
        "passed": overall_passed,
        "cases": summaries,
        "plots": [str(path.relative_to(repo_root)) if path.is_relative_to(repo_root) else str(path) for path in plot_paths],
        "notes": [
            "These are collision-only complete PICOS++ runs with field solve, RF, pair source, and axial motion disabled.",
            "Conservation checks here mean fixed super-particle population/weights and net charge in closed periodic cases.",
            "Thermal/equilibration tolerances are intentionally broad because these are small stochastic PIC runs.",
        ],
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0 if overall_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
