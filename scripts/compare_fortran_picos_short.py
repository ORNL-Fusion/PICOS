#!/usr/bin/env python3
"""Set up and compare a short LinearFokkerPlanck Fortran run with PICOS++.

The Fortran code is single-species. For this RF/ECH comparison it is configured
as kinetic electrons (`species_a = 1`). PICOS++ is run with ion and electron
species, and this script compares the PICOS++ electron species (`spp_2`) against
the Fortran electron output.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import struct
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import numpy as np


E_CHARGE = 1.602176634e-19
M_E = 9.1093837015e-31
C_LIGHT = 299792458.0


def parse_picos_input(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("//", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 2:
            values[parts[0]] = parts[1]
    return values


def rf_value(values: dict[str, str], species: str, key: str) -> str:
    species_key = f"RF_{species}_{key}"
    legacy_key = f"RF_{key}"
    if species_key in values:
        return values[species_key]
    return values[legacy_key]


def read_fortran_record_float64(path: Path) -> np.ndarray:
    payload_size = path.stat().st_size
    with path.open("rb") as handle:
        header = handle.read(4)
        if len(header) != 4:
            raise ValueError(f"{path} is too small to be a Fortran record")

        candidates = []
        for endian in ("<", ">"):
            nbytes = struct.unpack(f"{endian}I", header)[0]
            if nbytes > 0 and nbytes + 8 <= payload_size and nbytes % 8 == 0:
                candidates.append((endian, nbytes))
        if not candidates:
            raise ValueError(f"{path} has an invalid Fortran record header")

        endian, nbytes = candidates[0]
        data = np.fromfile(handle, dtype=f"{endian}f8", count=nbytes // 8)
        trailer = handle.read(4)
        if len(trailer) == 4:
            trailing_nbytes = struct.unpack(f"{endian}I", trailer)[0]
            if trailing_nbytes != nbytes:
                raise ValueError(f"{path} has mismatched Fortran record markers")
        return data


def h5dump_text(path: Path, dataset: str, header_only: bool = False) -> str:
    cmd = ["h5dump"]
    if header_only:
        cmd.append("-H")
    cmd.extend(["-d", dataset, str(path)])
    return subprocess.check_output(cmd, text=True)


def read_hdf5_dataset(path: Path, dataset: str) -> np.ndarray:
    header = h5dump_text(path, dataset, header_only=True)
    shape_match = re.search(r"DATASPACE\s+SIMPLE\s+\{\s+\(\s*([0-9,\s]+)\)", header)
    if not shape_match:
        raise ValueError(f"Could not parse HDF5 shape for {dataset} in {path}")
    shape = tuple(int(item.strip()) for item in shape_match.group(1).split(",") if item.strip())

    if "H5T_IEEE_F64LE" in header:
        dtype = "<f8"
    elif "H5T_IEEE_F32LE" in header:
        dtype = "<f4"
    else:
        raise ValueError(f"Unsupported HDF5 datatype in {dataset} from {path}")

    with tempfile.NamedTemporaryFile(prefix="picos_h5_", suffix=".bin") as payload, tempfile.NamedTemporaryFile(prefix="picos_h5_", suffix=".ddl") as ddl:
        subprocess.check_call(["h5dump", "-d", dataset, "-o", payload.name, "-b", "LE", "-O", ddl.name, str(path)])
        data = np.fromfile(payload.name, dtype=dtype)
    return data.reshape(shape)


def picos_components(velocity: np.ndarray) -> np.ndarray:
    if velocity.ndim != 2:
        raise ValueError(f"Expected 2D V_p dataset, got shape {velocity.shape}")
    if velocity.shape[0] in (2, 3):
        return velocity
    if velocity.shape[1] in (2, 3):
        return velocity.T
    raise ValueError(f"Cannot infer velocity component axis from shape {velocity.shape}")


def kinetic_energy_eV(speed: np.ndarray, relativistic: bool) -> np.ndarray:
    if not relativistic:
        return 0.5 * M_E * speed * speed / E_CHARGE
    beta2 = np.square(speed / C_LIGHT)
    gamma = 1.0 / np.sqrt(1.0 - np.clip(beta2, 0.0, 1.0 - 1.0e-15))
    return (gamma - 1.0) * M_E * C_LIGHT * C_LIGHT / E_CHARGE


def stats(values: np.ndarray, te_eV: float) -> dict[str, float]:
    selected = values[np.isfinite(values)]
    if selected.size == 0:
        return {key: math.nan for key in ("count", "mean", "std", "p50", "p90", "p95", "p99", "max", "frac_gt_10Te")}
    return {
        "count": float(selected.size),
        "mean": float(np.mean(selected)),
        "std": float(np.std(selected)),
        "p50": float(np.percentile(selected, 50)),
        "p90": float(np.percentile(selected, 90)),
        "p95": float(np.percentile(selected, 95)),
        "p99": float(np.percentile(selected, 99)),
        "max": float(np.max(selected)),
        "frac_gt_10Te": float(np.mean(selected > 10.0 * te_eV)),
    }


def make_fortran_case(args: argparse.Namespace, picos_tag: str) -> Path:
    linear_root = args.linear_root.resolve()
    picos_input = args.picos_root / "picosFILES" / "inputFiles" / f"input_file_{picos_tag}.input"
    values = parse_picos_input(picos_input)

    cv_b = float(values["CV_B"])
    z_min = float(values["LX_min"])
    z_max = float(values["LX_max"])
    n_profile = int(values["IC_BX_NX"])
    nx = n_profile - 2
    dx = (z_max - z_min) / nx
    z = z_min - 0.5 * dx + dx * np.arange(n_profile)

    b_norm = np.loadtxt(args.picos_root / "picosFILES" / "inputFiles" / values["IC_BX_fileName"])
    bfield_name = f"{args.fortran_case_name}_Bfield.txt"
    bfield_path = linear_root / "BfieldData" / bfield_name
    np.savetxt(bfield_path, np.column_stack([z, cv_b * b_norm]), fmt="%.16e")

    fortran_input = linear_root / "InputFiles" / f"{args.fortran_case_name}.in"
    fortran_text = f"""&params_nml
params%fileDescriptor = '{args.fortran_descriptor}',
params%repoDir       = "{linear_root}",
params%BFieldFileDir = "/BfieldData",
params%BFieldFile    = "/{bfield_name}",
params%nz            = {n_profile},

params%NC = {args.particles},
params%NS = {args.steps},
params%dt = {args.physical_time / args.steps:.16e},
params%G  = 0.0,

params%jstart = {args.steps},
params%jend   = {args.steps},
params%jincr  = 1,

params%dtheta = 6.283185307179586,
params%r1     = 0.0,
params%r2     = {float(values["r2"]):.16e},
params%zmax   = {z_max:.16e},
params%zmin   = {z_min:.16e},
params%NZmesh = {nx},

params%iSave      = .true.,
params%iPush      = .true.,
params%iColl      = .false.,
params%iHeat      = .true.,
params%iPotential = .false.,
params%iDrag      = .false.,

params%Te0          = {float(values["CV_Te"]):.16e},
params%Ti0          = {float(values["CV_Te"]):.16e},
params%ne0          = {float(values["CV_ne"]):.16e},
params%Aion         = 1.0,
params%Zion         = 1.0,
params%species_a    = 1,
params%elevel       = 10,
params%CollOperType = 2,

params%BC_Type    = 3,
params%BC_zp_mean = 0.0,
params%BC_zp_std  = 0.3,
params%BC_Ep      = 0.0,
params%BC_Tp      = {float(values["CV_Te"]):.16e},
params%BC_xip     = 0.707,

params%IC_Type    = 1,
params%IC_zp_mean = 0.0,
params%IC_zp_std  = 0.3,
params%IC_Ep      = 0.0,
params%IC_Tp      = {float(values["CV_Te"]):.16e},
params%IC_xip     = 0.707,

params%f_RF       = {float(rf_value(values, "electron", "freq")):.16e},
params%zRes1      = {float(rf_value(values, "electron", "x1")):.16e},
params%zRes2      = {float(rf_value(values, "electron", "x2")):.16e},
params%kper       = {float(rf_value(values, "electron", "kper")):.16e},
params%kpar       = {float(rf_value(values, "electron", "kpar")):.16e},
params%Prf        = {args.rf_power:.16e},
params%n_harmonic = {int(rf_value(values, "electron", "n_harmonic"))},

params%s1   = 0.0,
params%s2   = 1.3,
params%s3   = 4.0,
params%phi1 = 0.0,
params%phi2 = 0.0,
params%phi3 = 0.0
/
"""
    fortran_input.write_text(fortran_text)
    return fortran_input


def replace_picos_key(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"^({re.escape(key)}\s+)\S+.*$", re.MULTILINE)
    if not pattern.search(text):
        raise KeyError(f"Missing PICOS input key {key}")
    return pattern.sub(rf"\g<1>{value}", text)


def species_rf_block(values: dict[str, str], one_file: str, simulation_time: float, rf_power: float) -> str:
    def value(key: str) -> str:
        return rf_value(values, "electron", key)

    return f"""// Ion RF operator:
// =============================================================================
RF_ion_Prf                      {rf_power:.16e}
RF_ion_n_harmonic               {value("n_harmonic")}
RF_ion_freq                     {value("freq")}
RF_ion_x1                       {value("x1")}
RF_ion_x2                       {value("x2")}
RF_ion_t_ON                     0.0
RF_ion_t_OFF                    {simulation_time:.16e}
RF_ion_kpar                     {value("kpar")}
RF_ion_kper                     {value("kper")}
RF_ion_handedness               {value("handedness")}
RF_ion_EfieldMode               0
RF_ion_EfieldAmplitude          {value("EfieldAmplitude")}
RF_ion_maxEnergyGainFraction    {value("maxEnergyGainFraction")}
RF_ion_maxParticleEnergy        {value("maxParticleEnergy")}
RF_ion_maxVelocityFractionC     {value("maxVelocityFractionC")}
RF_ion_Prf_fileName             {one_file}
RF_ion_Prf_NS                   {values["IC_Te_NX"]}

// Electron RF/ECH operator:
// =============================================================================
RF_electron_Prf                      {rf_power:.16e}
RF_electron_n_harmonic               {value("n_harmonic")}
RF_electron_freq                     {value("freq")}
RF_electron_x1                       {value("x1")}
RF_electron_x2                       {value("x2")}
RF_electron_t_ON                     0.0
RF_electron_t_OFF                    {simulation_time:.16e}
RF_electron_kpar                     {value("kpar")}
RF_electron_kper                     {value("kper")}
RF_electron_handedness               {value("handedness")}
RF_electron_EfieldMode               0
RF_electron_EfieldAmplitude          {value("EfieldAmplitude")}
RF_electron_maxEnergyGainFraction    {value("maxEnergyGainFraction")}
RF_electron_maxParticleEnergy        {value("maxParticleEnergy")}
RF_electron_maxVelocityFractionC     {value("maxVelocityFractionC")}
RF_electron_Prf_fileName             {one_file}
RF_electron_Prf_NS                   {values["IC_Te_NX"]}

"""


def clone_local_picos_deck(args: argparse.Namespace, tag: str, relativistic: int) -> None:
    source_tag = "xray_case8_gc_rf_smoke"
    input_dir = args.picos_root / "picosFILES" / "inputFiles"
    source_input = input_dir / f"input_file_{source_tag}.input"
    source_ions = input_dir / f"ions_properties_{source_tag}.ion"
    if not source_input.is_file() or not source_ions.is_file():
        raise FileNotFoundError(f"Local fallback source deck {source_tag} is missing")

    source_values = parse_picos_input(source_input)
    source_simulation_time = float(source_values["simulationTime"])
    simulation_time = source_simulation_time * (args.physical_time / 2.0e-10)

    b_file = f"{tag}_B_norm.txt"
    one_file = f"{tag}_one.txt"
    shutil.copy2(input_dir / source_values["IC_BX_fileName"], input_dir / b_file)
    shutil.copy2(input_dir / source_values["IC_Te_fileName"], input_dir / one_file)

    input_text = source_input.read_text()
    input_text = input_text.replace(source_tag, tag)
    for key, value in {
        "quietStart": "0",
        "advanceParticleMethod": "1",
        "simulationTime": f"{simulation_time:.16e}",
        "SW_EfieldSolve": "0",
        "SW_fieldSolveModel": "0",
        "SW_BfieldSolve": "0",
        "SW_Collisions": "0",
        "SW_RFheating": "1",
        "SW_RFheatingIons": "0",
        "SW_RFheatingElectrons": "1",
        "SW_relativisticElectrons": str(relativistic),
        "IC_BX_fileName": b_file,
        "outputCadence": f"{simulation_time:.16e}",
    }.items():
        input_text = replace_picos_key(input_text, key, value)
    input_text = re.sub(
        r"// RF operator:\n// =+\n.*?\n(?=// Output variables:)",
        species_rf_block(source_values, one_file, simulation_time, args.rf_power),
        input_text,
        flags=re.S,
    )
    (input_dir / f"input_file_{tag}.input").write_text(input_text)

    ions_text = source_ions.read_text().replace(source_tag, tag)
    for key, value in {
        "NPC1": "2",
        "NPC2": "8",
        "pctSupPartOutput1": "1.0000000000000000e+02",
        "pctSupPartOutput2": "1.0000000000000000e+02",
    }.items():
        ions_text = replace_picos_key(ions_text, key, value)
    (input_dir / f"ions_properties_{tag}.ion").write_text(ions_text)


def generate_picos_deck(args: argparse.Namespace, tag: str, relativistic: int) -> None:
    local_source = args.picos_root / "picosFILES" / "inputFiles" / "input_file_xray_case8_gc_rf_smoke.input"
    if local_source.is_file():
        clone_local_picos_deck(args, tag, relativistic)
        return

    cmd = [
        "python3",
        "scripts/create_picos_xray_case.py",
        "--case",
        "Case8",
        "--tag",
        tag,
        "--physical-time",
        f"{args.physical_time:.16e}",
        "--nx",
        str(args.nx),
        "--ion-particles-per-cell",
        "2",
        "--electron-particles-per-cell",
        "8",
        "--output-saves",
        "1",
        "--collisions",
        "0",
        "--field-solve",
        "none",
        "--boundary-type",
        "3",
        "--advance-particle-method",
        "1",
        "--quiet-start",
        "0",
        "--rf-heat-ions",
        "0",
        "--rf-heat-electrons",
        "1",
        "--relativistic-electrons",
        str(relativistic),
        "--rf-efield-mode",
        "0",
        "--rf-power",
        f"{args.rf_power:.16e}",
        "--output-particle-percent",
        "100",
    ]
    subprocess.check_call(cmd, cwd=args.picos_root)


def read_fortran_output(args: argparse.Namespace) -> dict[str, np.ndarray]:
    output_dir = args.linear_root / "OutputFiles" / args.fortran_case_name / args.fortran_descriptor
    energy = read_fortran_record_float64(output_dir / "kep.out")
    z = read_fortran_record_float64(output_dir / "zp.out")
    pitch = read_fortran_record_float64(output_dir / "xip.out")
    energy_eV = energy.reshape((args.particles, -1), order="F")[:, -1]
    return {
        "energy_eV": energy_eV,
        "z_m": z.reshape((args.particles, -1), order="F")[:, -1],
        "pitch": pitch.reshape((args.particles, -1), order="F")[:, -1],
        "speed_over_c": np.sqrt(np.clip(2.0 * E_CHARGE * energy_eV / M_E, 0.0, None)) / C_LIGHT,
    }


def read_picos_output(args: argparse.Namespace, tag: str, relativistic: bool) -> dict[str, np.ndarray]:
    hdf_dir = args.picos_root / "picosFILES" / "outputFiles" / tag / "HDF5"
    velocities = []
    positions = []
    for path in sorted(hdf_dir.glob("PARTICLES_FILE_*.h5")):
        components = picos_components(read_hdf5_dataset(path, "/1/ions/spp_2/V_p"))
        velocities.append(components)
        x = read_hdf5_dataset(path, "/1/ions/spp_2/X_p").reshape(-1)
        positions.append(x)

    if not velocities:
        raise FileNotFoundError(f"No PICOS particle files found in {hdf_dir}")

    v = np.concatenate(velocities, axis=1)
    x = np.concatenate(positions)
    speed = np.sqrt(np.sum(v * v, axis=0))
    pitch = np.divide(v[0], speed, out=np.zeros_like(speed), where=speed > 0.0)
    return {
        "energy_eV": kinetic_energy_eV(speed, relativistic=relativistic),
        "z_m": x,
        "pitch": pitch,
        "speed_over_c": speed / C_LIGHT,
    }


def write_report(args: argparse.Namespace, rows: list[dict[str, Any]]) -> None:
    out_dir = args.picos_root / "validation" / "fortran_picos_compare"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "fortran_picos_short_compare.csv"
    md_path = out_dir / "fortran_picos_short_compare.md"

    fieldnames = list(rows[0].keys())
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# Fortran versus PICOS++ Short Comparison",
        "",
        f"Fortran case: `{args.linear_root / 'InputFiles' / (args.fortran_case_name + '.in')}`",
        f"PICOS++ tags: `{args.picos_tag_nonrel}`, `{args.picos_tag_rel}`",
        "",
        "| run | N | mean E [eV] | P95 [eV] | P99 [eV] | max E [eV] | mean z [m] | mean pitch | max v/c |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        max_vc = row.get("max_speed_over_c", math.nan)
        lines.append(
            f"| {row['run']} | {int(row['count'])} | {row['energy_mean_eV']:.4g} | "
            f"{row['energy_p95_eV']:.4g} | {row['energy_p99_eV']:.4g} | {row['energy_max_eV']:.4g} | "
            f"{row['z_mean_m']:.4g} | {row['pitch_mean']:.4g} | {max_vc:.4g} |"
        )

    by_run = {row["run"]: row for row in rows}
    if "fortran_nonrel" in by_run and "picos_nonrel" in by_run:
        delta = by_run["picos_nonrel"]["energy_mean_eV"] - by_run["fortran_nonrel"]["energy_mean_eV"]
        lines.extend(["", f"PICOS++ nonrel mean-energy delta from Fortran: `{delta:.6g} eV`."])
    if "fortran_nonrel" in by_run and "picos_rel" in by_run:
        delta = by_run["picos_rel"]["energy_mean_eV"] - by_run["fortran_nonrel"]["energy_mean_eV"]
        lines.append(f"PICOS++ relativistic mean-energy delta from Fortran: `{delta:.6g} eV`.")

    md_path.write_text("\n".join(lines) + "\n")
    print(md_path)
    print(csv_path)


def summarize_run(name: str, arrays: dict[str, np.ndarray], te_eV: float) -> dict[str, Any]:
    energy = stats(arrays["energy_eV"], te_eV)
    z = stats(arrays["z_m"], te_eV)
    pitch = stats(arrays["pitch"], te_eV)
    return {
        "run": name,
        "count": energy["count"],
        "energy_mean_eV": energy["mean"],
        "energy_std_eV": energy["std"],
        "energy_p50_eV": energy["p50"],
        "energy_p90_eV": energy["p90"],
        "energy_p95_eV": energy["p95"],
        "energy_p99_eV": energy["p99"],
        "energy_max_eV": energy["max"],
        "energy_frac_gt_10Te": energy["frac_gt_10Te"],
        "z_mean_m": z["mean"],
        "z_std_m": z["std"],
        "pitch_mean": pitch["mean"],
        "pitch_std": pitch["std"],
        "max_speed_over_c": float(np.max(arrays["speed_over_c"])) if "speed_over_c" in arrays else math.nan,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--picos-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--linear-root", type=Path, default=Path("/Users/78k/Desktop/picos_kinetic_electron_eval/LinearFokkerPlanck_Axisymmetric"))
    parser.add_argument("--fortran-case-name", default="xp_Case8PicosCompare")
    parser.add_argument("--fortran-descriptor", default="picos_compare_case8_power")
    parser.add_argument("--picos-tag-nonrel", default="xray_case8_fortran_compare_nonrel")
    parser.add_argument("--picos-tag-rel", default="xray_case8_fortran_compare_rel")
    parser.add_argument("--physical-time", type=float, default=2.0e-10)
    parser.add_argument("--steps", type=int, default=2141)
    parser.add_argument("--particles", type=int, default=640)
    parser.add_argument("--nx", type=int, default=80)
    parser.add_argument("--rf-power", type=float, default=3.0e5)
    parser.add_argument("--setup", action="store_true")
    parser.add_argument("--compare", action="store_true")
    args = parser.parse_args()

    args.picos_root = args.picos_root.resolve()
    args.linear_root = args.linear_root.resolve()

    if args.setup or not args.compare:
        generate_picos_deck(args, args.picos_tag_nonrel, relativistic=0)
        generate_picos_deck(args, args.picos_tag_rel, relativistic=1)
        make_fortran_case(args, args.picos_tag_nonrel)

    if args.compare:
        te_eV = float(parse_picos_input(args.picos_root / "picosFILES" / "inputFiles" / f"input_file_{args.picos_tag_nonrel}.input")["CV_Te"])
        rows = [
            summarize_run("fortran_nonrel", read_fortran_output(args), te_eV),
            summarize_run("picos_nonrel", read_picos_output(args, args.picos_tag_nonrel, relativistic=False), te_eV),
            summarize_run("picos_rel", read_picos_output(args, args.picos_tag_rel, relativistic=True), te_eV),
        ]
        write_report(args, rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
