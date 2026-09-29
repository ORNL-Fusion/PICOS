#!/usr/bin/env python3
"""Create matched local ECH-OFF/ECH-ON validation cases from one checkpoint."""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path


GYROPERIODS_20US = "2.0860734119736900e+02"
GYROPERIODS_5US = "5.2151835299342250e+01"
GYROPERIODS_0P5US = "5.2151835299342250e+00"
GYROPERIODS_1US = "1.0430367059868450e+01"
GYROPERIODS_0P1US = "1.0430367059868450e+00"


def replace(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"^({re.escape(key)}\s+).*$", re.MULTILINE)
    updated, count = pattern.subn(rf"\g<1>{value}", text)
    if count:
        return updated
    return text.rstrip() + f"\n{key:<42} {value}\n"


def configured(text: str, *, stage: str, restart_path: Path | None) -> str:
    common = {
        "Poisson_BCModel": "2",
        "SW_collisionConservationProjection": "1",
        "SW_kineticElectronQuasiNeutralProjection": "1",
        "kineticElectronQuasiNeutralIterations": "2",
        "kineticElectronQuasiNeutralRelaxation": "1.0",
        "SW_kineticElectronBackgroundHeating": "0",
        "kineticElectronBackgroundPower": "0.0",
        "RF_electron_resonanceMode": "0",
        "RF_electron_maxEnergyGainFraction": "0.0",
        "RF_electron_maxParticleEnergy": "0.0",
        "RF_electron_maxVelocityFractionC": "0.0",
    }
    for key, value in common.items():
        text = replace(text, key, value)

    if stage.startswith("profile_"):
        rf_on = stage == "profile_ech"
        values = {
            "simulationTime": GYROPERIODS_1US,
            "outputCadence": GYROPERIODS_0P1US,
            "SW_kineticElectronThermostat": "0",
            "kineticElectronThermostatRelaxation": "1.0",
            "SW_pairSource": "0",
            "SW_RFheating": "1" if rf_on else "0",
            "SW_RFheatingIons": "0",
            "SW_RFheatingElectrons": "1" if rf_on else "0",
            "restart_enabled": "0",
            "restart_path": "none",
            "restart_continueTime": "0",
            "RF_electron_Prf": "3.0000000000000000e+05",
            "RF_electron_freq": "7.0000000000000000e+10",
            "RF_electron_x1": "2.6000000000000001e+00",
            "RF_electron_x2": "2.7999999999999998e+00",
            "RF_electron_t_ON": "0.0",
            "RF_electron_t_OFF": GYROPERIODS_1US,
            "RF_electron_EfieldMode": "0",
        }
    elif stage == "precondition":
        values = {
            "simulationTime": GYROPERIODS_20US,
            "outputCadence": GYROPERIODS_5US,
            "SW_kineticElectronThermostat": "1",
            "kineticElectronThermostatRelaxation": "1.0",
            "SW_RFheating": "0",
            "SW_RFheatingIons": "0",
            "SW_RFheatingElectrons": "0",
            "restart_enabled": "0",
            "restart_path": "none",
            "restart_continueTime": "0",
        }
    else:
        if restart_path is None:
            raise ValueError("restart_path is required for comparison stages")
        rf_on = stage == "ech_on"
        values = {
            "simulationTime": GYROPERIODS_5US,
            "outputCadence": GYROPERIODS_0P5US,
            "SW_kineticElectronThermostat": "0",
            "kineticElectronThermostatRelaxation": "1.0",
            "SW_RFheating": "1" if rf_on else "0",
            "SW_RFheatingIons": "0",
            "SW_RFheatingElectrons": "1" if rf_on else "0",
            "restart_enabled": "1",
            "restart_path": str(restart_path),
            "restart_snapshot": "-1",
            # Reset the stage clock so the 0--5 us RF gate is reachable.
            "restart_continueTime": "0",
            "RF_electron_Prf": "3.0000000000000000e+05",
            "RF_electron_freq": "7.0000000000000000e+10",
            "RF_electron_x1": "2.6000000000000001e+00",
            "RF_electron_x2": "2.7999999999999998e+00",
            "RF_electron_t_ON": "0.0",
            "RF_electron_t_OFF": GYROPERIODS_5US,
            "RF_electron_EfieldMode": "0",
        }
        for key, value in values.items():
            text = replace(text, key, value)
        return text

    for key, value in values.items():
        text = replace(text, key, value)
    return text


def copy_case(source: Path, destination: Path, *, stage: str, restart_path: Path | None) -> None:
    input_source = source / "picosFILES" / "inputFiles"
    if not input_source.is_dir():
        input_source = source / "inputFiles"
    input_destination = destination / "picosFILES" / "inputFiles"
    input_destination.parent.mkdir(parents=True, exist_ok=False)
    shutil.copytree(input_source, input_destination)
    input_file = input_destination / "input_file.input"
    input_file.write_text(configured(input_file.read_text(), stage=stage, restart_path=restart_path))
    # Keep restart normalization and the otherwise-unused pre-restart marker
    # construction identical in all three stages.
    ion_file = input_destination / "ions_properties.ion"
    ion_text = replace(ion_file.read_text(), "NPC1", "41")
    ion_text = replace(ion_text, "NPC2", "165")
    if stage.startswith("profile_"):
        ion_text = replace(ion_text, "BC_type_1", "4")
        ion_text = replace(ion_text, "BC_type_2", "4")
        def one_file(pattern: str) -> str:
            matches = sorted(input_destination.glob(pattern))
            if len(matches) != 1:
                raise RuntimeError(f"Expected one {pattern} profile under {input_destination}, got {matches}")
            return matches[0].name

        b_name = one_file("*B_norm_from_hybrid.txt")
        ne_name = one_file("*ne_norm_from_hybrid.txt")
        te_name = one_file("*te_norm_from_hybrid.txt")
        ti_name = one_file("*ti_norm_from_hybrid.txt")
        input_text = replace(input_file.read_text(), "IC_BX_fileName", b_name)
        input_text = replace(input_text, "IC_Te_fileName", te_name)
        input_file.write_text(input_text)
        for suffix, name in (
            ("IC_Tper_fileName_1", ti_name),
            ("IC_Tpar_fileName_1", ti_name),
            ("IC_densityFraction_fileName_1", ne_name),
            ("IC_Tper_fileName_2", te_name),
            ("IC_Tpar_fileName_2", te_name),
            ("IC_densityFraction_fileName_2", ne_name),
        ):
            ion_text = replace(ion_text, suffix, name)
    ion_file.write_text(ion_text)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-root",
        type=Path,
        default=root / "validation" / "local_uncapped_steady_ech_validation",
    )
    args = parser.parse_args()
    out_root = args.out_root.resolve()
    if out_root.exists():
        raise FileExistsError(f"Refusing to overwrite {out_root}")

    steady_source = root / "validation/local_laptop_fixed_mpex_source_10us_sheathboth.3YgOHQ"
    ech_source = root / "validation/PICOS_NERSC_MPEX_kinetic_1ms_then_ech100us/left_resonance_70ghz"
    profile_source = root / (
        "validation/local_main_vs_fully_kinetic_profile_compare_10us_20260914/"
        "fully_kinetic_periodic_control_10us"
    )
    precondition = out_root / "precondition_20us"
    restart_path = precondition / "picosFILES/outputFiles/HDF5"
    copy_case(steady_source, precondition, stage="precondition", restart_path=None)
    copy_case(ech_source, out_root / "steady_hold_5us", stage="steady_hold", restart_path=restart_path)
    copy_case(ech_source, out_root / "ech_on_5us", stage="ech_on", restart_path=restart_path)
    copy_case(profile_source, out_root / "profile_control_1us", stage="profile_control", restart_path=None)
    copy_case(profile_source, out_root / "profile_ech_1us", stage="profile_ech", restart_path=None)

    print(out_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
