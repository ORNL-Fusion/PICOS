#!/usr/bin/env python3
"""Create ECH-OFF/ON cases initialized from the NF-validated hybrid steady state."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

from create_local_uncapped_steady_ech_validation import copy_case, replace


def read_h5(path: Path, dataset: str) -> np.ndarray:
    output = subprocess.run(
        ["h5dump", "-y", "-w", "0", "-d", dataset, str(path)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    match = re.search(r"DATA\s*\{(.*?)\}\s*\}\s*\}\s*$", output, re.S)
    if match is None:
        raise RuntimeError(f"Cannot parse {dataset} from {path}")
    return np.fromstring(match.group(1).replace(",", " "), sep=" ")


def hybrid_average(hybrid: Path, dataset: str, first: int, last: int) -> np.ndarray:
    particle_file = hybrid / "PARTICLES_FILE_0.h5"
    return np.mean(
        [read_h5(particle_file, f"/{step}/ions/spp_1/{dataset}") for step in range(first, last + 1)],
        axis=0,
    )


def configure_profiles(case: Path, hybrid: Path, first: int, last: int) -> None:
    inputs = case / "picosFILES/inputFiles"
    input_file = inputs / "input_file.input"
    ion_file = inputs / "ions_properties.ion"

    # quietStart previously selected a uniform initializer and silently ignored all
    # imported profiles.  The nonuniform initializer is mandatory for this gate.
    text = replace(input_file.read_text(), "quietStart", "0")
    input_file.write_text(text)
    ions = replace(ion_file.read_text(), "NPC1", "82")
    ions = replace(ions, "NPC2", "330")
    ion_file.write_text(ions)

    source_x = read_h5(hybrid / "main.h5", "/geometry/x_m")
    target_x = np.linspace(-2.0, 8.0, 200)
    density = hybrid_average(hybrid, "n_m", first, last)
    tpar = hybrid_average(hybrid, "Tpar_m", first, last)
    tper = hybrid_average(hybrid, "Tper_m", first, last)

    ne_file = next(inputs.glob("*ne_norm_from_hybrid.txt"))
    te_file = next(inputs.glob("*te_norm_from_hybrid.txt"))
    old_ti_file = next(inputs.glob("*ti_norm_from_hybrid.txt"))
    tpar_file = inputs / "NF_validated_Tpar_norm_from_hybrid.txt"
    tper_file = inputs / "NF_validated_Tper_norm_from_hybrid.txt"
    np.savetxt(ne_file, np.interp(target_x, source_x, density) / 1.0e19, fmt="%.16e")
    np.savetxt(te_file, np.ones_like(target_x), fmt="%.16e")
    np.savetxt(tpar_file, np.interp(target_x, source_x, tpar) / 15.0, fmt="%.16e")
    np.savetxt(tper_file, np.interp(target_x, source_x, tper) / 15.0, fmt="%.16e")
    old_ti_file.unlink()

    ions = ion_file.read_text()
    ions = replace(ions, "IC_Tpar_fileName_1", tpar_file.name)
    ions = replace(ions, "IC_Tper_fileName_1", tper_file.name)
    ion_file.write_text(ions)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument(
        "--hybrid",
        type=Path,
        default=Path("/Users/78k/Desktop/PICOS_fully_kinetic_matched_validation_20260925/hybrid"),
    )
    parser.add_argument("--first", type=int, default=10)
    parser.add_argument("--last", type=int, default=22)
    args = parser.parse_args()
    out = args.out_root.resolve()
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")

    source = root / (
        "validation/local_main_vs_fully_kinetic_profile_compare_10us_20260914/"
        "fully_kinetic_periodic_control_10us"
    )
    for name, stage in (("profile_control_1us", "profile_control"), ("profile_ech_1us", "profile_ech")):
        case = out / name
        copy_case(source, case, stage=stage, restart_path=None)
        configure_profiles(case, args.hybrid.resolve(), args.first, args.last)
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
