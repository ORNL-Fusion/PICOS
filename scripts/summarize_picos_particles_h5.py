#!/usr/bin/env python3
"""Summarize PICOS particle velocities from an HDF5 output file.

This script intentionally uses the HDF5 command-line tools instead of h5py so
it works on the current Mac environment without installing Python packages.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np


E_CHARGE = 1.602176634e-19
M_E = 9.1093837015e-31
C_LIGHT = 299792458.0


def h5dump_text(path: Path, dataset: str, header_only: bool = False) -> str:
    cmd = ["h5dump"]
    if header_only:
        cmd.append("-H")
    cmd.extend(["-d", dataset, str(path)])
    return subprocess.check_output(cmd, text=True)


def read_dataset_raw(path: Path, dataset: str) -> np.ndarray:
    header = h5dump_text(path, dataset, header_only=True)
    shape_match = re.search(r"DATASPACE\s+SIMPLE\s+\{\s+\(\s*([0-9,\s]+)\)", header)
    if not shape_match:
        raise ValueError(f"Could not parse HDF5 shape for {dataset}")
    shape = tuple(int(item.strip()) for item in shape_match.group(1).split(",") if item.strip())

    if "H5T_IEEE_F64LE" in header:
        dtype = "<f8"
    elif "H5T_IEEE_F32LE" in header:
        dtype = "<f4"
    else:
        raise ValueError(f"Unsupported HDF5 datatype in {dataset}")

    with tempfile.NamedTemporaryFile(prefix="picos_h5_", suffix=".bin") as payload, tempfile.NamedTemporaryFile(prefix="picos_h5_", suffix=".ddl") as ddl:
        subprocess.check_call(["h5dump", "-d", dataset, "-o", payload.name, "-b", "LE", "-O", ddl.name, str(path)])
        data = np.fromfile(payload.name, dtype=dtype)
    return data.reshape(shape)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("particles_file", type=Path)
    parser.add_argument("--iteration", default="1")
    parser.add_argument("--species", default="spp_2")
    parser.add_argument("--mass-kg", type=float, default=M_E)
    parser.add_argument("--relativistic", action="store_true", help="Use KE=(gamma-1)mc^2 instead of nonrelativistic 0.5*m*v^2.")
    args = parser.parse_args()

    dataset = f"/{args.iteration}/ions/{args.species}/V_p"
    velocity = read_dataset_raw(args.particles_file, dataset)
    if velocity.ndim != 2:
        raise ValueError(f"Expected 2D V_p dataset, got shape {velocity.shape}")

    if velocity.shape[0] in (2, 3):
        components = velocity
    elif velocity.shape[1] in (2, 3):
        components = velocity.T
    else:
        raise ValueError(f"Cannot infer velocity component axis from shape {velocity.shape}")

    speed = np.sqrt(np.sum(components * components, axis=0))
    beta2 = np.square(speed / C_LIGHT)
    superluminal_particles = int(np.count_nonzero(beta2 >= 1.0))
    gamma = 1.0 / np.sqrt(1.0 - np.clip(beta2, 0.0, 1.0 - 1.0e-15))
    nonrel_energy_eV = 0.5 * args.mass_kg * speed * speed / E_CHARGE
    rel_energy_eV = (gamma - 1.0) * args.mass_kg * C_LIGHT * C_LIGHT / E_CHARGE
    energy_eV = rel_energy_eV if args.relativistic else nonrel_energy_eV
    summary = {
        "dataset": dataset,
        "shape": list(velocity.shape),
        "energy_model": "relativistic" if args.relativistic else "nonrelativistic",
        "finite": bool(np.isfinite(velocity).all()),
        "particles": int(speed.size),
        "superluminal_particles": superluminal_particles,
        "max_speed_m_s": float(np.max(speed)),
        "max_speed_over_c": float(np.max(speed) / C_LIGHT),
        "mean_energy_eV": float(np.mean(energy_eV)),
        "p50_energy_eV": float(np.percentile(energy_eV, 50)),
        "p95_energy_eV": float(np.percentile(energy_eV, 95)),
        "p99_energy_eV": float(np.percentile(energy_eV, 99)),
        "max_energy_eV": float(np.max(energy_eV)),
        "mean_energy_nonrel_eV": float(np.mean(nonrel_energy_eV)),
        "mean_energy_rel_eV": float(np.mean(rel_energy_eV)),
        "max_energy_nonrel_eV": float(np.max(nonrel_energy_eV)),
        "max_energy_rel_eV": float(np.max(rel_energy_eV)),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
