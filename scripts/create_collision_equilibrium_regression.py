#!/usr/bin/env python3
"""Create homogeneous, field-free kinetic electron-ion collision regressions."""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import numpy as np


def replace_key(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(key)}\s+).*$")
    if not pattern.search(text):
        return text + f"\n{key:<36} {value}\n"
    return pattern.sub(rf"\g<1>{value}", text)


def configure_case(case: Path, duration_us: float, cadence_us: float,
                   self_collisions: int, cross_collisions: int,
                   global_projection: int) -> None:
    inputs = case / "picosFILES" / "inputFiles"
    output = case / "picosFILES" / "outputFiles"
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)

    input_path = inputs / "input_file.input"
    text = input_path.read_text()
    gyroperiods_per_us = 1.0430367059868450e1
    updates = {
        "simulationTime": f"{duration_us * gyroperiods_per_us:.16e}",
        "outputCadence": f"{cadence_us * gyroperiods_per_us:.16e}",
        # Keep the kinetic field configuration internally valid.  With the
        # particle push disabled, the solved field cannot change velocities.
        "SW_EfieldSolve": "1",
        "SW_Collisions": "1",
        "SW_RFheating": "0",
        "SW_RFheatingIons": "0",
        "SW_RFheatingElectrons": "0",
        "SW_pairSource": "0",
        "SW_advancePos": "0",
        "SW_Bohm": "0",
        "SW_kineticElectronQuasiNeutralProjection": "0",
        "SW_kineticElectronThermostat": "0",
        "SW_kineticElectronBackgroundHeating": "0",
        "SW_collisionSelfSpecies": str(self_collisions),
        "SW_collisionCrossSpecies": str(cross_collisions),
        "SW_collisionSelfConservation": "1",
        "SW_collisionConservationProjection": str(global_projection),
        "IC_uniformBfield": "1",
        "IC_BX": "1.0",
        "IC_ne": "1.0000000000000000e+20",
        "CV_ne": "1.0000000000000000e+20",
        "CV_B": "1.0",
        "CV_Te": "1.5000000000000000e+01",
        "CV_Tpar": "1.5000000000000000e+01",
        "CV_Tper": "1.5000000000000000e+01",
    }
    for key, value in updates.items():
        text = replace_key(text, key, value)
    input_path.write_text(text)

    ion_path = inputs / "ions_properties.ion"
    ions = ion_path.read_text()
    for species in (1, 2):
        ions = replace_key(ions, f"IC_Tpar_{species}", "1.5000000000000000e+01")
        ions = replace_key(ions, f"IC_Tper_{species}", "1.5000000000000000e+01")
        ions = replace_key(ions, f"IC_Upar_{species}", "0.0")
    ion_path.write_text(ions)

    # Every referenced initial-condition shape is made spatially uniform.  The
    # test then exercises collisions only: no field push, sources, streaming,
    # or boundary replacement can change either species' energy.
    for profile in inputs.glob("*.txt"):
        if any(token in profile.name for token in ("ne_norm", "te_norm", "Tpar_norm", "Tper_norm")):
            np.savetxt(profile, np.ones(200), fmt="%.16e")
        elif "Upar_norm" in profile.name:
            np.savetxt(profile, np.zeros(200), fmt="%.16e")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--duration-us", type=float, default=2.0)
    parser.add_argument("--cadence-us", type=float, default=0.1)
    args = parser.parse_args()
    if args.out_root.exists():
        raise FileExistsError(f"Refusing to overwrite {args.out_root}")

    cases = {
        "none": (0, 0, 0),
        "self_only": (1, 0, 0),
        "cross_only_unprojected": (0, 1, 0),
        "cross_only_projected": (0, 1, 1),
        "full_projected": (1, 1, 1),
    }
    for name, switches in cases.items():
        case = args.out_root / name
        shutil.copytree(args.base, case)
        configure_case(case, args.duration_us, args.cadence_us, *switches)
    print(args.out_root.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
