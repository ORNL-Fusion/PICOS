#!/usr/bin/env python3
"""Summarize archived 2020 Proto-MPEX X-ray/ECH benchmark outputs.

The archived LinearFokkerPlanck_Axisymmetric outputs are Fortran unformatted
records: a uint32 byte count, float64 payload, and a trailing uint32 byte count.
This script extracts particle histories, computes final-energy summary
statistics, locates the cyclotron resonance field crossing, and writes simple
CSV/PNG comparison artifacts for PICOS validation.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import struct
from pathlib import Path
from typing import Any

import numpy as np


E_CHARGE = 1.602176634e-19
M_E = 9.1093837015e-31


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
        key, value = match.group(1), match.group(2).strip().rstrip(",")
        values[key] = parse_value(value)

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


def particle_history(values: np.ndarray, n_particles: int, n_saves: int) -> np.ndarray:
    expected = n_particles * n_saves
    if values.size != expected:
        inferred = values.size // max(n_particles, 1)
        if inferred * n_particles != values.size:
            raise ValueError(
                f"Cannot reshape {values.size} values into particle history with {n_particles} particles"
            )
        n_saves = inferred
    return values.reshape((n_particles, n_saves), order="F")


def find_bfield_file(case_dir: Path, xray_root: Path, metadata: dict[str, Any]) -> Path | None:
    name = str(metadata.get("BFieldFile", "")).strip("/")
    candidates = [
        case_dir / name,
        xray_root / "BfieldData" / name,
        case_dir / Path(name).name,
        xray_root / "BfieldData" / Path(name).name,
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def resonance_locations(bfield_path: Path | None, rf_frequency_hz: float, harmonic: int) -> tuple[float, str]:
    if harmonic == 0:
        return math.nan, ""
    b_res = 2.0 * math.pi * rf_frequency_hz * M_E / (E_CHARGE * abs(harmonic))
    if bfield_path is None:
        return b_res, ""

    data = np.loadtxt(bfield_path)
    z = data[:, 0]
    b = data[:, 1]
    order = np.argsort(z)
    z = z[order]
    b = b[order]

    roots: list[float] = []
    diff = b - b_res
    for i in range(len(diff) - 1):
        if diff[i] == 0:
            roots.append(float(z[i]))
        elif diff[i] * diff[i + 1] < 0:
            f = abs(diff[i]) / (abs(diff[i]) + abs(diff[i + 1]))
            roots.append(float(z[i] + f * (z[i + 1] - z[i])))
    return b_res, ";".join(f"{root:.6g}" for root in roots)


def safe_stats(values: np.ndarray, clip_eV: float, te_eV: float) -> dict[str, float]:
    finite = np.isfinite(values)
    physical = finite & (values >= 0.0)
    clipped = physical & (values <= clip_eV)
    selected = values[clipped]

    stats: dict[str, float] = {
        "valid": float(np.count_nonzero(physical)),
        "nonfinite": float(values.size - np.count_nonzero(finite)),
        "above_clip": float(np.count_nonzero(physical & (values > clip_eV))),
    }
    if selected.size == 0:
        for key in ("mean", "p50", "p90", "p95", "p99", "max", "frac_gt_10Te", "frac_gt_100eV", "frac_gt_500eV", "frac_gt_1000eV"):
            stats[key] = math.nan
        return stats

    stats["mean"] = float(np.mean(selected))
    stats["p50"] = float(np.percentile(selected, 50))
    stats["p90"] = float(np.percentile(selected, 90))
    stats["p95"] = float(np.percentile(selected, 95))
    stats["p99"] = float(np.percentile(selected, 99))
    stats["max"] = float(np.max(selected))
    stats["frac_gt_10Te"] = float(np.mean(selected > 10.0 * te_eV))
    stats["frac_gt_100eV"] = float(np.mean(selected > 100.0))
    stats["frac_gt_500eV"] = float(np.mean(selected > 500.0))
    stats["frac_gt_1000eV"] = float(np.mean(selected > 1000.0))
    return stats


def summarize_case(case_dir: Path, xray_root: Path, clip_eV: float) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    metadata = parse_case_input(case_dir / "xp_Xray.in")
    n_particles = int(metadata["Nparts"])
    times = read_fortran_record_float64(case_dir / "tp.out")
    n_saves = int(times.size)

    kinetic_energy = particle_history(read_fortran_record_float64(case_dir / "kep.out"), n_particles, n_saves)
    z_position = particle_history(read_fortran_record_float64(case_dir / "zp.out"), n_particles, n_saves)
    pitch = particle_history(read_fortran_record_float64(case_dir / "xip.out"), n_particles, n_saves)

    final_energy = kinetic_energy[:, -1]
    final_z = z_position[:, -1]
    final_pitch = pitch[:, -1]
    te_eV = float(metadata.get("Te0", math.nan))
    stats = safe_stats(final_energy, clip_eV, te_eV)

    bfield_path = find_bfield_file(case_dir, xray_root, metadata)
    b_res, z_res = resonance_locations(
        bfield_path,
        float(metadata.get("f_RF", math.nan)),
        int(metadata.get("n_harmonic", 0)),
    )

    row: dict[str, Any] = {
        "case": case_dir.name,
        "heated": bool(metadata.get("iHeat", False)),
        "potential": bool(metadata.get("iPotential", False)),
        "n_particles": n_particles,
        "n_steps": int(metadata.get("Nsteps", -1)),
        "saved_frames": n_saves,
        "dt_s": float(metadata.get("dt", math.nan)),
        "t_final_s": float(times[-1]) if times.size else math.nan,
        "ne0_m3": float(metadata.get("ne0", math.nan)),
        "Te0_eV": te_eV,
        "Ti0_eV": float(metadata.get("Ti0", math.nan)),
        "rf_frequency_Hz": float(metadata.get("f_RF", math.nan)),
        "harmonic": int(metadata.get("n_harmonic", 0)),
        "Ew_Vm": float(metadata.get("Ew", math.nan)),
        "kpar_m-1": float(metadata.get("kpar", math.nan)),
        "kper_m-1": float(metadata.get("kper", math.nan)),
        "Bres_T": b_res,
        "z_resonance_m": z_res,
        "Bfield_file": str(bfield_path) if bfield_path else "",
    }
    row.update({f"final_energy_{key}_eV": value for key, value in stats.items()})

    for name in ("pcount1", "pcount2", "pcount3", "pcount4", "ecount1", "ecount2", "ecount3", "ecount4"):
        path = case_dir / f"{name}.out"
        if path.is_file():
            arr = read_fortran_record_float64(path)
            row[f"{name}_sum"] = float(np.nansum(arr))
            row[f"{name}_last"] = float(arr[-1]) if arr.size else math.nan

    arrays = {
        "time_s": times,
        "final_energy_eV": final_energy,
        "final_z_m": final_z,
        "final_pitch": final_pitch,
    }
    return row, arrays


def write_plots(out_dir: Path, rows: list[dict[str, Any]], arrays_by_case: dict[str, dict[str, np.ndarray]], clip_eV: float) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
    except Exception as exc:  # pragma: no cover - optional plotting dependency
        print(f"Skipping plots because matplotlib is unavailable: {exc}")
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    for row in rows:
        case = row["case"]
        energy = arrays_by_case[case]["final_energy_eV"]
        selected = energy[np.isfinite(energy) & (energy >= 0.0) & (energy <= clip_eV)]
        if selected.size:
            upper = min(clip_eV, max(50.0, float(np.percentile(selected, 99.7))))
            bins = np.linspace(0.0, upper, 120)
            ax.hist(selected, bins=bins, density=True, histtype="step", linewidth=1.8, label=case)
    ax.set_xlabel("Final electron energy [eV]")
    ax.set_ylabel("PDF")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "ech_xray_final_energy_hist.png", dpi=200)
    plt.close(fig)

    heated_cases = [row["case"] for row in rows if row["heated"]]
    if heated_cases:
        case = heated_cases[0]
        energy = arrays_by_case[case]["final_energy_eV"]
        z = arrays_by_case[case]["final_z_m"]
        selected = np.isfinite(energy) & np.isfinite(z) & (energy >= 0.0) & (energy <= clip_eV)
        if np.count_nonzero(selected) > 0:
            fig, ax = plt.subplots(figsize=(8, 5))
            upper = min(clip_eV, max(50.0, float(np.percentile(energy[selected], 99.5))))
            image = ax.hist2d(z[selected], energy[selected], bins=(120, 120), range=[[np.min(z[selected]), np.max(z[selected])], [0.0, upper]], norm=LogNorm())
            ax.set_xlabel("z [m]")
            ax.set_ylabel("Final electron energy [eV]")
            ax.set_title(f"{case}: final EEDF versus z")
            fig.colorbar(image[3], ax=ax, label="particles/bin")
            fig.tight_layout()
            fig.savefig(out_dir / "ech_xray_energy_vs_z.png", dpi=200)
            plt.close(fig)


def write_markdown(out_dir: Path, rows: list[dict[str, Any]]) -> None:
    lines = [
        "# ECH X-ray Benchmark Summary",
        "",
        "Primary comparison target: heated Case8 versus no-heating Case11 from the archived 2020 Proto-MPEX X-ray study.",
        "",
        "| Case | Heated | Mean E [eV] | P95 [eV] | P99 [eV] | >10Te | >500 eV | Bres [T] | z_res [m] |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            "| {case} | {heated} | {mean:.3g} | {p95:.3g} | {p99:.3g} | {f10:.3g} | {f500:.3g} | {bres:.4g} | {zres} |".format(
                case=row["case"],
                heated=int(bool(row["heated"])),
                mean=row["final_energy_mean_eV"],
                p95=row["final_energy_p95_eV"],
                p99=row["final_energy_p99_eV"],
                f10=row["final_energy_frac_gt_10Te_eV"],
                f500=row["final_energy_frac_gt_500eV_eV"],
                bres=row["Bres_T"],
                zres=row["z_resonance_m"],
            )
        )

    by_case = {row["case"]: row for row in rows}
    if "Case8" in by_case and "Case11" in by_case:
        heated = by_case["Case8"]
        cold = by_case["Case11"]
        lines.extend(
            [
                "",
                "## Heated versus No-heating Delta",
                "",
                f"Case8 mean final energy / Case11 mean final energy: {heated['final_energy_mean_eV'] / cold['final_energy_mean_eV']:.3g}",
                f"Case8 P95 - Case11 P95: {heated['final_energy_p95_eV'] - cold['final_energy_p95_eV']:.3g} eV",
                f"Case8 fast fraction (>10Te) - Case11 fast fraction: {heated['final_energy_frac_gt_10Te_eV'] - cold['final_energy_frac_gt_10Te_eV']:.3g}",
            ]
        )

    (out_dir / "ech_xray_summary.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xray-root", type=Path, default=Path("/Users/78k/Desktop/2020_07_13_XrayStudy"))
    parser.add_argument("--cases", nargs="+", default=["Case8", "Case11", "Case15"])
    parser.add_argument("--out-dir", type=Path, default=Path("validation/xray_study"))
    parser.add_argument("--energy-clip-eV", type=float, default=1.0e6)
    args = parser.parse_args()

    xray_root = args.xray_root.expanduser().resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    arrays_by_case: dict[str, dict[str, np.ndarray]] = {}
    for case in args.cases:
        case_dir = xray_root / "OutputFiles" / "xp_Xray" / case
        row, arrays = summarize_case(case_dir, xray_root, args.energy_clip_eV)
        rows.append(row)
        arrays_by_case[row["case"]] = arrays

    csv_path = out_dir / "ech_xray_case_summary.csv"
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    write_markdown(out_dir, rows)
    write_plots(out_dir, rows, arrays_by_case, args.energy_clip_eV)

    print(f"Wrote {csv_path}")
    print(f"Wrote {out_dir / 'ech_xray_summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
