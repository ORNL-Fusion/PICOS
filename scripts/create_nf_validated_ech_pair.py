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


def configure_profiles(case: Path, hybrid: Path, first: int, last: int,
                       ion_npc: int, electron_npc: int,
                       cross_collisions: bool) -> None:
    inputs = case / "picosFILES/inputFiles"
    input_file = inputs / "input_file.input"
    ion_file = inputs / "ions_properties.ion"

    # quietStart previously selected a uniform initializer and silently ignored all
    # imported profiles.  The nonuniform initializer is mandatory for this gate.
    text = replace(input_file.read_text(), "quietStart", "0")
    # Open kinetic runs need charge-balanced end losses.  A fixed 3 Te cutoff
    # selectively removes the far electron tail without guaranteeing that the
    # transmitted electron current matches the ion current; over many transit
    # times that cools the EEDF and eventually opens charge-density holes.
    text = replace(text, "Poisson_sheathCurrentBalance", "1")
    text = replace(text, "SW_Bohm", "1")
    text = replace(text, "Bohm_type", "2")
    text = replace(text, "Bohm_edgeCells", "4")
    text = replace(text, "Bohm_t_ON", "0.0")
    text = replace(text, "Bohm_gamma_i", "3.0")
    text = replace(text, "SW_collisionSelfSpecies", "1")
    text = replace(text, "SW_collisionSelfConservation", "1")
    text = replace(text, "SW_collisionCrossSpecies", "1" if cross_collisions else "0")
    input_file.write_text(text)
    ions = replace(ion_file.read_text(), "NPC1", str(ion_npc))
    ions = replace(ions, "NPC2", str(electron_npc))
    ion_file.write_text(ions)

    source_x = read_h5(hybrid / "main.h5", "/geometry/x_m")
    target_x = np.linspace(-2.0, 8.0, 200)
    density = hybrid_average(hybrid, "n_m", first, last)
    tpar = hybrid_average(hybrid, "Tpar_m", first, last)
    tper = hybrid_average(hybrid, "Tper_m", first, last)
    upar = hybrid_average(hybrid, "u_m/x", first, last)

    ne_file = next(inputs.glob("*ne_norm_from_hybrid.txt"))
    te_file = next(inputs.glob("*te_norm_from_hybrid.txt"))
    old_ti_file = next(inputs.glob("*ti_norm_from_hybrid.txt"))
    tpar_file = inputs / "NF_validated_Tpar_norm_from_hybrid.txt"
    tper_file = inputs / "NF_validated_Tper_norm_from_hybrid.txt"
    upar_file = inputs / "NF_validated_Upar_norm_from_hybrid.txt"
    np.savetxt(ne_file, np.interp(target_x, source_x, density) / 1.0e19, fmt="%.16e")
    np.savetxt(te_file, np.ones_like(target_x), fmt="%.16e")
    np.savetxt(tpar_file, np.interp(target_x, source_x, tpar) / 15.0, fmt="%.16e")
    np.savetxt(tper_file, np.interp(target_x, source_x, tper) / 15.0, fmt="%.16e")
    np.savetxt(upar_file, np.interp(target_x, source_x, upar) / 1.0e5, fmt="%.16e")
    old_ti_file.unlink()

    ions = ion_file.read_text()
    ions = replace(ions, "IC_Tpar_fileName_1", tpar_file.name)
    ions = replace(ions, "IC_Tper_fileName_1", tper_file.name)
    for species in (1, 2):
        ions = replace(ions, f"IC_Upar_{species}", "1.0000000000000000e+05")
        ions = replace(ions, f"IC_Upar_fileName_{species}", upar_file.name)
        ions = replace(ions, f"IC_Upar_NX_{species}", "200")
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
    parser.add_argument("--duration-us", type=float, default=1.0)
    parser.add_argument("--cadence-us", type=float, default=0.1)
    # Boundary-current balance is quantized by the ion macro-particle charge.
    # 82 ions/cell made one ion loss worth roughly 16 electron markers and
    # produced bursty energy losses; 330/cell reduces that ratio to about four.
    parser.add_argument("--ion-npc", type=int, default=330)
    # Roughly 1,300 electron markers/cell are needed here to resolve the small
    # tail that crosses the downstream ambipolar-potential barrier.  The prior
    # 330/cell local shortcut depleted that tail and produced a false field
    # runaway after a few microseconds.
    parser.add_argument("--electron-npc", type=int, default=1320)
    parser.add_argument("--cross-collisions", action="store_true",
                        help="Enable kinetic electron-ion collisions (off by default for an RF-only isolation pair)")
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
        configure_profiles(case, args.hybrid.resolve(), args.first, args.last,
                           args.ion_npc, args.electron_npc, args.cross_collisions)
        input_file = case / "picosFILES/inputFiles/input_file.input"
        text = input_file.read_text()
        gyroperiods_per_us = 1.0430367059868450e1
        duration = args.duration_us*gyroperiods_per_us
        cadence = args.cadence_us*gyroperiods_per_us
        text = replace(text, "simulationTime", f"{duration:.16e}")
        text = replace(text, "outputCadence", f"{cadence:.16e}")
        text = replace(text, "RF_electron_t_OFF", f"{duration:.16e}")
        input_file.write_text(text)
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
