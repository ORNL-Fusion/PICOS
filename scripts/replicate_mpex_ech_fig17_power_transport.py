#!/usr/bin/env python3
"""Set up and analyze Fig. 17-style MPEX ECH target power transport.

This is not the short ECH operator benchmark.  Fig. 17-style transport requires
energy reaching the dump/target boundaries, so the comparison quantity is the
boundary energy flux from Fortran `ecount1/ecount2` and PICOS++ HDF5
`/boundary/E1` and `/boundary/E2`.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from compare_fortran_picos_short import (
    ensure_plotting,
    generate_picos_deck,
    make_fortran_case,
    parse_picos_input,
    read_fortran_record_float64,
    read_hdf5_scalar_optional,
    replace_picos_key,
    replace_or_insert_picos_key,
)


ELECTRON_CFL_DT_S = 9.342275e-14
ELECTRON_THERMAL_SPEED_16EV = 2.37e6


@dataclass
class CaseSpec:
    power_w: float
    fortran_case: str
    fortran_descriptor: str
    picos_tag_nonrel: str
    picos_tag_rel: str


def power_label(power_w: float) -> str:
    return f"{int(round(power_w / 1000.0)):03d}kW"


def parse_power_list(text: str) -> list[float]:
    values: list[float] = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        values.append(float(item) * 1000.0)
    if not values:
        raise ValueError("At least one RF power must be supplied")
    return values


def case_specs(args: argparse.Namespace) -> list[CaseSpec]:
    specs: list[CaseSpec] = []
    for power_w in parse_power_list(args.rf_powers_kw):
        label = power_label(power_w)
        specs.append(
            CaseSpec(
                power_w=power_w,
                fortran_case=f"{args.case_prefix}_{label}",
                fortran_descriptor=f"{args.descriptor_prefix}_{label}",
                picos_tag_nonrel=f"{args.case_prefix.lower()}_{label.lower()}_nonrel",
                picos_tag_rel=f"{args.case_prefix.lower()}_{label.lower()}_rel",
            )
        )
    return specs


def patch_fortran_transport_case(path: Path, args: argparse.Namespace) -> None:
    text = path.read_text()
    replacements = {
        "params%BC_Type": str(args.boundary_type),
        "params%G": f"{args.source_rate:.16e}",
        "params%BC_zp_mean": f"{args.source_z_mean:.16e}",
        "params%BC_zp_std": f"{args.source_z_sigma:.16e}",
        "params%IC_zp_mean": f"{args.source_z_mean:.16e}",
        "params%IC_zp_std": f"{args.source_z_sigma:.16e}",
    }
    for key, value in replacements.items():
        text = replace_fortran_key(text, key, value)
    path.write_text(text)


def replace_fortran_key(text: str, key: str, value: str) -> str:
    import re

    pattern = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*)[^,\n]+(,?)", re.MULTILINE)
    if not pattern.search(text):
        raise KeyError(f"Missing Fortran input key {key}")
    return pattern.sub(rf"\g<1>{value}\2", text, count=1)


def patch_picos_transport_case(tag: str, args: argparse.Namespace) -> None:
    input_dir = args.picos_root / "picosFILES" / "inputFiles"
    input_path = input_dir / f"input_file_{tag}.input"
    ions_path = input_dir / f"ions_properties_{tag}.ion"

    input_text = input_path.read_text()
    header_marker = "// =============================================================================\nmpisForFields"
    header = f"""// PICOS input generated for MPEX ECH Fig. 17-style power transport
// Source case: 2020 X-ray/Proto-MPEX ECH Case8 profile and RF parameters.
// Model notes:
// - Particles: 1D-2V guiding-center D+ ions plus kinetic electrons.
// - Field: electrostatic field solve is disabled for this transport scaffold.
// - Heating: electron RF/ECH is enabled; ion RF heating is disabled.
// - Boundaries: open/reinjecting source boundaries are used so target/dump power can be measured.
// - Relativity: {'enabled' if tag.endswith('_rel') else 'disabled'} for the kinetic-electron particle energy and RF paths.
// Field solver model legend:
// - 0 = Ohm-law hybrid, 1 = charge-density Poisson, 2 = reformulated Poisson stress-moment E-field update.
"""
    if header_marker in input_text:
        input_text = header + input_text.split(header_marker, 1)[1].join(["// =============================================================================\nmpisForFields", ""])
    input_text = replace_or_insert_picos_key(
        input_text,
        "SW_relativisticRFOperator",
        "1" if tag.endswith("_rel") else "0",
        "SW_relativisticElectrons",
    )
    input_text = replace_or_insert_picos_key(input_text, "IC_randomSeed", "271828", "IC_velocityDistributionModel")
    input_text = replace_or_insert_picos_key(input_text, "CollOperType", "2", "SW_Collisions")
    input_text = replace_or_insert_picos_key(input_text, "SW_collisionConservationProjection", "0", "CollOperType")
    input_text = replace_or_insert_picos_key(input_text, "collisionRandomSeed", "314159", "SW_collisionConservationProjection")
    input_text = replace_or_insert_picos_key(input_text, "SW_pairSource", "0", "SW_RFheatingElectrons")
    input_path.write_text(input_text)

    values = parse_picos_input(input_path)
    one_file = values["IC_Te_fileName"]
    ion_ppc = max(1, int(round(args.particles / args.nx / 4.0)))
    electron_ppc = max(1, int(round(args.particles / args.nx)))

    ions_text = ions_path.read_text()
    ion_header_marker = "// =============================================================================\n// Species 1"
    ion_header = f"""// PICOS species generated for MPEX ECH Fig. 17-style power transport
// Used by input_file_{tag}.input.
// Species are 1D-2V guiding-center D+ ions and kinetic electrons.
// Boundary type 1 is used by default to record and reinject particles at dump/target boundaries.
"""
    if ion_header_marker in ions_text:
        ions_text = ion_header + ions_text.split(ion_header_marker, 1)[1].join(["// =============================================================================\n// Species 1", ""])
    for key, value in {
        "pairSource_ionSpecies": "1",
        "pairSource_electronSpecies": "2",
        "pairSource_rate": "0.0",
        "pairSource_mean_x": f"{args.source_z_mean:.16e}",
        "pairSource_sigma_x": f"{args.source_z_sigma:.16e}",
        "pairSource_Ti_birth": "1.6000000000000000e+01",
        "pairSource_Te_birth": "1.6000000000000000e+01",
        "pairSource_Ei_birth": "0",
        "pairSource_Ee_birth": "0",
        "pairSource_eta_i": "0",
        "pairSource_eta_e": "0",
        "pairSource_positionMode": "0",
        "pairSource_fileName": one_file,
        "pairSource_NS": values["IC_Te_NX"],
        "pairSource_maxParticleWeight": "1000",
    }.items():
        input_text = replace_or_insert_picos_key(input_text, key, value, "IC_Te_fileName")
    input_path.write_text(input_text)

    for key, value in {
        "NPC1": str(ion_ppc),
        "NPC2": str(electron_ppc),
        "pctSupPartOutput1": "1.0000000000000000e+02",
        "pctSupPartOutput2": "1.0000000000000000e+02",
        "BC_type_1": str(args.boundary_type),
        "BC_type_2": str(args.boundary_type),
        "BC_mean_x_1": f"{args.source_z_mean:.16e}",
        "BC_mean_x_2": f"{args.source_z_mean:.16e}",
        "BC_sigma_x_1": f"{args.source_z_sigma:.16e}",
        "BC_sigma_x_2": f"{args.source_z_sigma:.16e}",
        "BC_G_1": f"{args.ion_source_rate:.16e}",
        "BC_G_2": f"{args.source_rate:.16e}",
        "BC_G_fileName_1": one_file,
        "BC_G_fileName_2": one_file,
        "BC_G_NS_1": values["IC_Te_NX"],
        "BC_G_NS_2": values["IC_Te_NX"],
    }.items():
        ions_text = replace_picos_key(ions_text, key, value)
    ions_path.write_text(ions_text)


def setup_cases(args: argparse.Namespace) -> list[CaseSpec]:
    specs = case_specs(args)
    for spec in specs:
        local_args = argparse.Namespace(**vars(args))
        local_args.rf_power = spec.power_w
        local_args.fortran_case_name = spec.fortran_case
        local_args.fortran_descriptor = spec.fortran_descriptor
        local_args.picos_tag_nonrel = spec.picos_tag_nonrel
        local_args.picos_tag_rel = spec.picos_tag_rel

        generate_picos_deck(local_args, spec.picos_tag_nonrel, relativistic=0)
        generate_picos_deck(local_args, spec.picos_tag_rel, relativistic=1)
        patch_picos_transport_case(spec.picos_tag_nonrel, args)
        patch_picos_transport_case(spec.picos_tag_rel, args)

        fortran_input = make_fortran_case(local_args, spec.picos_tag_nonrel)
        patch_fortran_transport_case(fortran_input, args)
    return specs


def fortran_output_dir(args: argparse.Namespace, spec: CaseSpec) -> Path:
    return args.linear_root / "OutputFiles" / spec.fortran_case / spec.fortran_descriptor


def picos_hdf_dir(args: argparse.Namespace, tag: str) -> Path:
    return args.picos_root / "picosFILES" / "outputFiles" / tag / "HDF5"


def latest_particle_h5(hdf_dir: Path) -> Path | None:
    files = sorted(hdf_dir.glob("PARTICLES_FILE_*.h5"))
    if not files:
        return None
    return files[0]


def finite_last(values: np.ndarray) -> float:
    values = np.ravel(np.asarray(values, dtype=float))
    values = values[np.isfinite(values)]
    if values.size == 0:
        return math.nan
    return float(values[-1])


def finite_window_mean(values: np.ndarray, fraction: float) -> float:
    values = np.ravel(np.asarray(values, dtype=float))
    values = values[np.isfinite(values)]
    if values.size == 0:
        return math.nan
    fraction = min(max(float(fraction), 0.0), 1.0)
    count = max(1, int(math.ceil(values.size * fraction)))
    return float(np.mean(values[-count:]))


def read_fortran_boundary(args: argparse.Namespace, spec: CaseSpec) -> dict[str, float]:
    output_dir = fortran_output_dir(args, spec)
    row: dict[str, float] = {
        "dump_power_kw": math.nan,
        "target_power_kw": math.nan,
        "rf_absorbed_kw": math.nan,
        "dump_particle_rate_s": math.nan,
        "target_particle_rate_s": math.nan,
    }
    if not output_dir.is_dir():
        return row

    files = {
        "dump_power_kw": "ecount1.out",
        "target_power_kw": "ecount2.out",
        "rf_absorbed_kw": "ecount3.out",
        "dump_particle_rate_s": "pcount1.out",
        "target_particle_rate_s": "pcount2.out",
    }
    for key, name in files.items():
        path = output_dir / name
        if not path.is_file():
            continue
        try:
            value = finite_window_mean(read_fortran_record_float64(path), args.average_tail_fraction)
        except Exception:
            continue
        row[key] = value / 1000.0 if key.endswith("_kw") else value
    return row


def read_picos_boundary(args: argparse.Namespace, tag: str) -> dict[str, float]:
    hdf_dir = picos_hdf_dir(args, tag)
    h5 = latest_particle_h5(hdf_dir)
    row: dict[str, float] = {
        "dump_power_kw": math.nan,
        "target_power_kw": math.nan,
        "rf_absorbed_kw": math.nan,
        "dump_particle_rate_s": math.nan,
        "target_particle_rate_s": math.nan,
    }
    if h5 is None:
        return row

    row["dump_power_kw"] = read_hdf5_scalar_optional(h5, "/1/boundary/E1") / 1000.0
    row["target_power_kw"] = read_hdf5_scalar_optional(h5, "/1/boundary/E2") / 1000.0
    row["rf_absorbed_kw"] = read_hdf5_scalar_optional(h5, "/1/rf/electron/E3") / 1000.0
    row["dump_particle_rate_s"] = read_hdf5_scalar_optional(h5, "/1/boundary/N1")
    row["target_particle_rate_s"] = read_hdf5_scalar_optional(h5, "/1/boundary/N2")
    return row


def collect_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for spec in case_specs(args):
        base = {
            "applied_rf_power_kw": spec.power_w / 1000.0,
            "fortran_case": spec.fortran_case,
            "fortran_descriptor": spec.fortran_descriptor,
        }
        for run, data in (
            ("fortran_nonrel", read_fortran_boundary(args, spec)),
            ("picos_nonrel", read_picos_boundary(args, spec.picos_tag_nonrel)),
            ("picos_rel", read_picos_boundary(args, spec.picos_tag_rel)),
        ):
            target = data["target_power_kw"]
            rows.append(
                {
                    **base,
                    "run": run,
                    "picos_tag": "" if run == "fortran_nonrel" else (spec.picos_tag_nonrel if run == "picos_nonrel" else spec.picos_tag_rel),
                    "dump_power_kw": data["dump_power_kw"],
                    "target_power_kw": target,
                    "rf_absorbed_kw": data["rf_absorbed_kw"],
                    "target_efficiency": target / (spec.power_w / 1000.0) if spec.power_w > 0.0 and math.isfinite(target) else math.nan,
                    "dump_particle_rate_s": data["dump_particle_rate_s"],
                    "target_particle_rate_s": data["target_particle_rate_s"],
                    "has_output": all(math.isfinite(data[key]) for key in ("dump_power_kw", "target_power_kw")),
                }
            )
    baseline: dict[str, dict[str, float]] = {}
    for row in rows:
        if float(row["applied_rf_power_kw"]) == 0.0 and row["has_output"]:
            baseline[row["run"]] = {
                "target_power_kw": float(row["target_power_kw"]),
                "dump_power_kw": float(row["dump_power_kw"]),
            }
    for row in rows:
        base = baseline.get(row["run"])
        if base is None or not row["has_output"]:
            row["target_increment_kw"] = math.nan
            row["dump_increment_kw"] = math.nan
        else:
            row["target_increment_kw"] = float(row["target_power_kw"]) - base["target_power_kw"]
            row["dump_increment_kw"] = float(row["dump_power_kw"]) - base["dump_power_kw"]
    return rows


def write_scan_csv(args: argparse.Namespace, rows: list[dict[str, Any]]) -> Path:
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "mpex_ech_fig17_power_transport_scan.csv"
    fieldnames = [
        "applied_rf_power_kw",
        "run",
        "target_power_kw",
        "target_increment_kw",
        "dump_power_kw",
        "dump_increment_kw",
        "rf_absorbed_kw",
        "target_efficiency",
        "dump_particle_rate_s",
        "target_particle_rate_s",
        "has_output",
        "fortran_case",
        "fortran_descriptor",
        "picos_tag",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def write_commands(args: argparse.Namespace) -> Path:
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "run_commands.sh"
    picos_files = args.picos_root / "picosFILES"
    xpicos = args.picos_root / "build" / "picosFILES" / "src" / "xpicos"
    linfp = args.linear_root / "src" / "linFP"
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "RUN_FORTRAN=${RUN_FORTRAN:-1}",
        "RUN_PICOS_TRANSPORT=${RUN_PICOS_TRANSPORT:-0}",
        "",
        "if [[ \"${RUN_FORTRAN}\" == \"1\" ]]; then",
        "# Fortran reference runs",
    ]
    for spec in case_specs(args):
        lines.extend(
            [
                f"cd {args.linear_root}",
                "OMP_NUM_THREADS=1 OMP_PROC_BIND=false \\",
                f"REPO_DIR={args.linear_root} \\",
                f"INPUT_FILE={spec.fortran_case}.in \\",
                f"INPUT_FILE_DIR={args.linear_root / 'InputFiles' / (spec.fortran_case + '.in')} \\",
                f"{linfp}",
                "",
            ]
        )
    lines.extend(
        [
            "fi",
            "",
            "if [[ \"${RUN_PICOS_TRANSPORT}\" == \"1\" ]]; then",
            "# PICOS++ kinetic-electron transport runs.",
            "# These are intentionally opt-in because the electron-CFL step makes",
            "# physical MPEX transport-time runs very expensive.",
        ]
    )
    for spec in case_specs(args):
        for tag in (spec.picos_tag_nonrel, spec.picos_tag_rel):
            lines.extend(
                [
                    f"cd {picos_files}",
                    "HDF5_USE_FILE_LOCKING=FALSE OMP_NUM_THREADS=1 \\",
                    f"mpirun -np {args.mpi_ranks} {xpicos} 1-D outputFiles {tag}",
                    "",
                ]
            )
    lines.extend(
        [
            "else",
            "echo \"Skipping full PICOS++ transport runs. Set RUN_PICOS_TRANSPORT=1 to launch them.\"",
            "fi",
            "",
            "# Replot after the runs finish",
            f"cd {args.picos_root}",
            "PYTHONPYCACHEPREFIX=/private/tmp/picos_pycache MPLCONFIGDIR=/private/tmp/picos_mpl \\",
            "python3 scripts/replicate_mpex_ech_fig17_power_transport.py --plot",
            "",
        ]
    )
    path.write_text("\n".join(lines))
    path.chmod(0o755)
    return path


def plot_rows(args: argparse.Namespace, rows: list[dict[str, Any]]) -> list[Path]:
    plt = ensure_plotting()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []

    powers = np.linspace(0.0, max(parse_power_list(args.rf_powers_kw)) / 1000.0, 200)
    fig, ax = plt.subplots(figsize=(13, 7), constrained_layout=True)
    ax.plot(powers, args.reference_efficiency * powers, color="0.25", ls="--", lw=2.0, label=f"{100.0 * args.reference_efficiency:.0f}% reference trend")
    ax.scatter([70.0], [14.0], color="0.1", marker="*", s=130, label="reported ~14 kW at 70 kW")

    styles = {
        "fortran_nonrel": ("k", "o", "Fortran nonrel"),
        "picos_nonrel": ("tab:blue", "s", "PICOS++ nonrel"),
        "picos_rel": ("tab:red", "^", "PICOS++ rel"),
    }
    for run, (color, marker, label) in styles.items():
        selected = [row for row in rows if row["run"] == run and row["has_output"]]
        if not selected:
            continue
        x = [float(row["applied_rf_power_kw"]) for row in selected]
        y = [float(row["target_power_kw"]) for row in selected]
        ax.plot(x, y, color=color, marker=marker, lw=1.8, ms=6.0, label=label)

    ax.set_xlabel("applied 28 GHz ECH power [kW]")
    ax.set_ylabel("power to target boundary [kW]")
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.28)
    ax.legend()
    path = out_dir / "mpex_ech_fig17_target_power.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(13, 7), constrained_layout=True)
    ax.plot(powers, args.reference_efficiency * powers, color="0.25", ls="--", lw=2.0, label=f"{100.0 * args.reference_efficiency:.0f}% reference trend")
    ax.scatter([70.0], [14.0], color="0.1", marker="*", s=130, label="reported ~14 kW at 70 kW")
    plotted_increment = False
    for run, (color, marker, label) in styles.items():
        selected = [row for row in rows if row["run"] == run and row["has_output"] and math.isfinite(row["target_increment_kw"])]
        if not selected:
            continue
        plotted_increment = True
        x = [float(row["applied_rf_power_kw"]) for row in selected]
        y = [float(row["target_increment_kw"]) for row in selected]
        ax.plot(x, y, color=color, marker=marker, lw=1.8, ms=6.0, label=label)
    ax.axhline(0.0, color="0.55", lw=1.0)
    ax.set_xlabel("applied 28 GHz ECH power [kW]")
    ax.set_ylabel("target power increment from 0 kW case [kW]")
    ax.set_xlim(left=0.0)
    if plotted_increment:
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(min(ymin, -5.0), max(ymax, 16.0))
    else:
        ax.set_ylim(-1.0, 16.0)
    ax.grid(True, alpha=0.28)
    ax.legend()
    path = out_dir / "mpex_ech_fig17_target_increment.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(13, 7), constrained_layout=True)
    plotted_boundary = False
    for run, (color, marker, label) in styles.items():
        selected = [row for row in rows if row["run"] == run and row["has_output"]]
        if not selected:
            continue
        plotted_boundary = True
        x = [float(row["applied_rf_power_kw"]) for row in selected]
        target = [float(row["target_power_kw"]) for row in selected]
        dump = [float(row["dump_power_kw"]) for row in selected]
        ax.plot(x, target, color=color, marker=marker, lw=1.8, ms=6.0, label=f"{label} target")
        ax.plot(x, dump, color=color, marker=marker, lw=1.2, ms=5.0, ls=":", label=f"{label} dump")
    ax.set_xlabel("applied 28 GHz ECH power [kW]")
    ax.set_ylabel("boundary energy flux [kW]")
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.28)
    if plotted_boundary:
        ax.legend(ncols=2, fontsize=24)
    else:
        ax.text(
            0.5,
            0.5,
            "No completed transport outputs found yet",
            ha="center",
            va="center",
            transform=ax.transAxes,
            color="0.35",
        )
    path = out_dir / "mpex_ech_fig17_target_dump_power.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    return paths


def write_report(args: argparse.Namespace, rows: list[dict[str, Any]], csv_path: Path, commands_path: Path, plots: list[Path]) -> Path:
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    physical_time = args.physical_time
    steps = physical_time / ELECTRON_CFL_DT_S
    transit_time = (args.target_z - args.source_z_mean) / ELECTRON_THERMAL_SPEED_16EV
    transit_steps = transit_time / ELECTRON_CFL_DT_S
    path = out_dir / "README.md"
    lines = [
        "# MPEX ECH Fig. 17 Power Transport Workflow",
        "",
        "This directory contains the Fig. 17-style target/dump power transport workflow.",
        "It is separate from `validation/mpex_ech_fig17_proxy`, which validates the local ECH operator using particle energy distributions.",
        "",
        "## What is compared",
        "",
        "- Fortran dump power: `ecount1.out`.",
        f"- Fortran target power: `ecount2.out`, averaged over the final `{100.0 * args.average_tail_fraction:.0f}%` of saved time steps.",
        "- PICOS++ dump power: HDF5 dataset `/1/boundary/E1`.",
        "- PICOS++ target power: HDF5 dataset `/1/boundary/E2`.",
        "- PICOS++ absorbed RF power diagnostic: HDF5 dataset `/1/rf/electron/E3`.",
        "- The `target_increment_kw` diagnostic subtracts the matching 0 kW case and is the quantity to compare against additional ECH transport to the target.",
        "",
        "The right boundary (`LX_max`) is treated as the target and the left boundary (`LX_min`) is treated as the dump, matching the existing PICOS boundary comments and the Fortran target-at-positive-z convention.",
        "",
        "## Important limitation",
        "",
        "The requested local PDF could not be read in this Codex session because macOS denied file access with `Operation not permitted`.",
        "The plot therefore uses the accessible published trend information as a reference: target load of order 14 kW at 70 kW applied 28 GHz power and about 18% transport efficiency.",
        "",
        "A full kinetic-electron PICOS++ reproduction of target power is much more expensive than the short ECH operator test:",
        "",
        f"- Current requested physical time for this deck family: `{physical_time:.6e} s`, about `{steps:.3e}` electron-CFL steps.",
        f"- Estimated source-to-target transit time from `z={args.source_z_mean:.3g} m` to `z={args.target_z:.3g} m` at 16 eV thermal speed: `{transit_time:.3e} s`, about `{transit_steps:.3e}` electron-CFL steps.",
        "- The short `2e-10 s` validation run cannot reproduce Fig. 17 target power because particles have not had transport time to reach the target/dump boundaries.",
        "",
        "## Generated files",
        "",
        f"- Scan table: `{csv_path.name}`.",
        f"- Run commands: `{commands_path.name}`.",
    ]
    for plot in plots:
        lines.append(f"- Plot: `{plot.name}`.")
    completed = [row for row in rows if row["has_output"]]
    if completed:
        lines.extend(
            [
                "",
                "## Current Extracted Scan",
                "",
                "| run | applied RF [kW] | target total [kW] | target increment [kW] | dump total [kW] | RF absorbed [kW] |",
                "|---|---:|---:|---:|---:|---:|",
            ]
        )
        for row in completed:
            lines.append(
                f"| {row['run']} | {float(row['applied_rf_power_kw']):.3g} | "
                f"{float(row['target_power_kw']):.4g} | {float(row['target_increment_kw']):.4g} | "
                f"{float(row['dump_power_kw']):.4g} | {float(row['rf_absorbed_kw']):.4g} |"
            )
        lines.extend(
            [
                "",
                "The present Fortran scan is not yet a validated Fig. 17 reproduction: the raw thermal boundary flux is much larger than the kW-scale reference, and the baseline-subtracted target increment is noisy rather than monotonic.",
            ]
        )
    lines.extend(
        [
            "",
            "## Run sequence",
            "",
            "```bash",
            "cd /Users/78k/Desktop/picos_kinetic_electron_eval/PICOS",
            "PYTHONPYCACHEPREFIX=/private/tmp/picos_pycache \\",
            "python3 scripts/replicate_mpex_ech_fig17_power_transport.py --setup --plot",
            "",
            "# Runs the Fortran scan and skips full PICOS++ transport by default.",
            "validation/mpex_ech_fig17_power_transport/run_commands.sh",
            "",
            "# Full PICOS++ transport is opt-in because it is tens of millions of steps.",
            "RUN_PICOS_TRANSPORT=1 validation/mpex_ech_fig17_power_transport/run_commands.sh",
            "",
            "PYTHONPYCACHEPREFIX=/private/tmp/picos_pycache MPLCONFIGDIR=/private/tmp/picos_mpl \\",
            "python3 scripts/replicate_mpex_ech_fig17_power_transport.py --plot",
            "```",
            "",
            "For production comparison, run the Fortran transport cases first, then choose a PICOS++ strategy that is computationally feasible while preserving the electron time-scale requirement.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")
    return path


def run_command(cmd: list[str], cwd: Path, env: dict[str, str] | None = None) -> None:
    subprocess.check_call(cmd, cwd=cwd, env=env)


def run_fortran_cases(args: argparse.Namespace) -> None:
    import os

    env_base = os.environ.copy()
    env_base["OMP_NUM_THREADS"] = "1"
    env_base["OMP_PROC_BIND"] = "false"
    env_base["REPO_DIR"] = str(args.linear_root)
    for spec in case_specs(args):
        env = env_base.copy()
        env["INPUT_FILE"] = f"{spec.fortran_case}.in"
        env["INPUT_FILE_DIR"] = str(args.linear_root / "InputFiles" / f"{spec.fortran_case}.in")
        run_command([str(args.linear_root / "src" / "linFP")], cwd=args.linear_root, env=env)


def run_picos_cases(args: argparse.Namespace) -> None:
    import os

    env = os.environ.copy()
    env["HDF5_USE_FILE_LOCKING"] = "FALSE"
    env["OMP_NUM_THREADS"] = "1"
    xpicos = args.picos_root / "build" / "picosFILES" / "src" / "xpicos"
    picos_files = args.picos_root / "picosFILES"
    for spec in case_specs(args):
        for tag in (spec.picos_tag_nonrel, spec.picos_tag_rel):
            run_command(["mpirun", "-np", str(args.mpi_ranks), str(xpicos), "1-D", "outputFiles", tag], cwd=picos_files, env=env)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--picos-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--linear-root", type=Path, default=Path("/Users/78k/Desktop/picos_kinetic_electron_eval/LinearFokkerPlanck_Axisymmetric"))
    parser.add_argument("--out-dir", type=Path, default=Path("validation/mpex_ech_fig17_power_transport"))
    parser.add_argument("--case-prefix", default="xp_MPEX_ECH_Fig17Transport")
    parser.add_argument("--descriptor-prefix", default="fig17_power_transport")
    parser.add_argument("--rf-powers-kw", default="0,10,30,50,70")
    parser.add_argument("--physical-time", type=float, default=5.0e-6)
    parser.add_argument("--steps", type=int, default=10000)
    parser.add_argument("--particles", type=int, default=6400)
    parser.add_argument("--nx", type=int, default=80)
    parser.add_argument("--source-rate", type=float, default=1.0e19)
    parser.add_argument("--ion-source-rate", type=float, default=1.0e19)
    parser.add_argument("--source-z-mean", type=float, default=0.0)
    parser.add_argument("--source-z-sigma", type=float, default=0.3)
    parser.add_argument("--target-z", type=float, default=4.3)
    parser.add_argument("--boundary-type", type=int, default=1)
    parser.add_argument("--mpi-ranks", type=int, default=4)
    parser.add_argument("--average-tail-fraction", type=float, default=0.5)
    parser.add_argument("--reference-efficiency", type=float, default=0.18)
    parser.add_argument("--setup", action="store_true")
    parser.add_argument("--run-fortran", action="store_true")
    parser.add_argument("--run-picos", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    args.picos_root = args.picos_root.resolve()
    args.linear_root = args.linear_root.resolve()
    if not args.out_dir.is_absolute():
        args.out_dir = args.picos_root / args.out_dir
    args.out_dir = args.out_dir.resolve()

    if args.setup:
        setup_cases(args)
    commands_path = write_commands(args)

    if args.run_fortran:
        run_fortran_cases(args)
    if args.run_picos:
        run_picos_cases(args)

    rows = collect_rows(args)
    csv_path = write_scan_csv(args, rows)
    plots: list[Path] = []
    if args.plot or not (args.setup or args.run_fortran or args.run_picos):
        plots = plot_rows(args, rows)
    report_path = write_report(args, rows, csv_path, commands_path, plots)

    print(report_path)
    print(csv_path)
    for plot in plots:
        print(plot)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
