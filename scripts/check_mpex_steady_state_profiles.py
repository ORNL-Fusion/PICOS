#!/usr/bin/env python3
"""Check whether a PICOS++ MPEX steady-state run has settled.

The script reads PICOS HDF5 output, plots plasma profiles versus time, and
compares late-time profile windows.  It is intended for the RF-off MPEX steady
stage before restarting into ECH.
"""

from __future__ import annotations

import argparse
import csv
import math
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SPECIES_LABEL = {
    "spp_1": "D+",
    "spp_2": "electron",
}


@dataclass
class SteadyCase:
    name: str
    hdf5_dir: Path
    z_m: np.ndarray
    times_s: np.ndarray
    species: list[str]
    mesh: dict[str, dict[str, np.ndarray]]
    fields: dict[str, np.ndarray]


def numeric_keys(handle: h5py.File) -> list[str]:
    return sorted([key for key in handle.keys() if key.isdigit()], key=lambda item: int(item))


def find_hdf5_dirs(path: Path, extract_dir: Path | None = None) -> list[Path]:
    path = path.expanduser()
    if path.is_file() and path.suffixes[-2:] == [".tar", ".gz"]:
        with tarfile.open(path, "r:gz") as tar:
            hdf5_main = [member for member in tar.getmembers() if member.name.endswith("/HDF5/main.h5")]
            if not hdf5_main:
                return []
            if extract_dir is None:
                extract_dir = Path("/private/tmp") / f"picos_steady_extract_{path.stem.replace('.tar', '')}"
            extract_dir.mkdir(parents=True, exist_ok=True)
            tar.extractall(extract_dir)
            return find_hdf5_dirs(extract_dir)

    if path.is_file() and path.name == "main.h5":
        return [path.parent.resolve()]
    if path.is_dir() and (path / "main.h5").is_file():
        return [path.resolve()]
    if path.is_dir():
        return sorted({candidate.parent.resolve() for candidate in path.rglob("HDF5/main.h5")})
    raise FileNotFoundError(path)


def read_scalar(handle: h5py.File, dataset: str, default: float = math.nan) -> float:
    if dataset not in handle:
        return default
    values = np.asarray(handle[dataset])
    return float(values.reshape(-1)[0]) if values.size else default


def dataset_1d(handle: h5py.File, dataset: str) -> np.ndarray | None:
    if dataset not in handle:
        return None
    return np.asarray(handle[dataset], dtype=float).reshape(-1)


def read_main_geometry(hdf5_dir: Path) -> np.ndarray:
    with h5py.File(hdf5_dir / "main.h5", "r") as handle:
        if "geometry/x_m" not in handle:
            raise KeyError(f"{hdf5_dir}/main.h5 has no geometry/x_m")
        return np.asarray(handle["geometry/x_m"], dtype=float).reshape(-1)


def read_species_names(hdf5_dir: Path, steps: list[str]) -> list[str]:
    p0 = hdf5_dir / "PARTICLES_FILE_0.h5"
    with h5py.File(p0, "r") as handle:
        root = handle[f"{steps[0]}/ions"]
        return sorted(root.keys(), key=lambda item: int(item.split("_")[-1]))


def read_case(hdf5_dir: Path) -> SteadyCase:
    p0 = hdf5_dir / "PARTICLES_FILE_0.h5"
    if not p0.is_file():
        raise FileNotFoundError(f"Expected particle-root output file: {p0}")

    with h5py.File(p0, "r") as handle:
        steps = numeric_keys(handle)
        if not steps:
            raise ValueError(f"No numeric snapshots in {p0}")
        times_s = np.asarray([read_scalar(handle, f"{step}/time") for step in steps], dtype=float)
        species = read_species_names(hdf5_dir, steps)
        mesh: dict[str, dict[str, np.ndarray]] = {}
        for spp in species:
            mesh[spp] = {}
            for var in ("n_m", "Tpar_m", "Tper_m"):
                cols = []
                for step in steps:
                    values = dataset_1d(handle, f"{step}/ions/{spp}/{var}")
                    if values is None:
                        cols = []
                        break
                    cols.append(values)
                if cols:
                    mesh[spp][var] = np.column_stack(cols)

            cols = []
            for step in steps:
                values = dataset_1d(handle, f"{step}/ions/{spp}/u_m/x")
                if values is None:
                    values = dataset_1d(handle, f"{step}/ions/{spp}/x")
                if values is None:
                    cols = []
                    break
                cols.append(values)
            if cols:
                mesh[spp]["u_m"] = np.column_stack(cols)

    fields: dict[str, np.ndarray] = {}
    f0 = hdf5_dir / "FIELDS_FILE_0.h5"
    if f0.is_file():
        with h5py.File(f0, "r") as handle:
            for var in ("BX_m", "EX_m", "Phi_m"):
                cols = []
                for step in steps:
                    values = dataset_1d(handle, f"{step}/fields/{var}/x")
                    if values is None:
                        cols = []
                        break
                    cols.append(values)
                if cols:
                    fields[var] = np.column_stack(cols)

    z_m = read_main_geometry(hdf5_dir)
    n = min([z_m.size] + [arr.shape[0] for per_spp in mesh.values() for arr in per_spp.values()] + [arr.shape[0] for arr in fields.values()])
    return SteadyCase(
        name=hdf5_dir.parent.name if hdf5_dir.parent.name != "outputFiles" else hdf5_dir.name,
        hdf5_dir=hdf5_dir,
        z_m=z_m[:n],
        times_s=times_s,
        species=species,
        mesh={spp: {var: arr[:n, :] for var, arr in values.items()} for spp, values in mesh.items()},
        fields={var: arr[:n, :] for var, arr in fields.items()},
    )


def smooth_z(values: np.ndarray, window: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if window <= 1:
        return values
    kernel = np.ones(window, dtype=float) / float(window)
    if values.ndim == 1:
        return np.convolve(values, kernel, mode="same")
    return np.apply_along_axis(lambda col: np.convolve(col, kernel, mode="same"), 0, values)


def relative_l2(a: np.ndarray, b: np.ndarray, floor: float = 1.0e-300) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if not np.any(mask):
        return math.nan
    denom = max(float(np.linalg.norm(b[mask])), floor)
    return float(np.linalg.norm(a[mask] - b[mask]) / denom)


def finite_mean_profile(values: np.ndarray, start: int, end: int) -> np.ndarray:
    return np.nanmean(values[:, start:end], axis=1)


def window_indices(n_times: int, window_fraction: float) -> tuple[slice, slice]:
    window = max(1, int(round(n_times * window_fraction)))
    window = min(window, max(1, n_times // 2))
    final = slice(n_times - window, n_times)
    previous = slice(n_times - 2 * window, n_times - window)
    return previous, final


def metric_rows(case: SteadyCase, window_fraction: float, tolerance: float) -> list[dict[str, object]]:
    if case.times_s.size < 4:
        return [
            {
                "case": case.name,
                "quantity": "all",
                "status": "not_enough_snapshots",
                "n_snapshots": case.times_s.size,
                "reason": "Need at least 4 output snapshots for a late-window steady-state check.",
            }
        ]

    prev_slice, final_slice = window_indices(case.times_s.size, window_fraction)
    rows: list[dict[str, object]] = []
    for spp in case.species:
        label = SPECIES_LABEL.get(spp, spp)
        for var in ("n_m", "Tpar_m", "Tper_m", "u_m"):
            if var not in case.mesh.get(spp, {}):
                continue
            values = case.mesh[spp][var]
            prev = finite_mean_profile(values, prev_slice.start, prev_slice.stop)
            final = finite_mean_profile(values, final_slice.start, final_slice.stop)
            last_change = relative_l2(values[:, -1], values[:, -2])
            window_change = relative_l2(prev, final)
            integral = float(np.trapz(final, case.z_m)) if final.size == case.z_m.size else math.nan
            peak_idx = int(np.nanargmax(final)) if np.any(np.isfinite(final)) else 0
            rows.append(
                {
                    "case": case.name,
                    "species": label,
                    "variable": var,
                    "status": "steady_like" if window_change <= tolerance and last_change <= tolerance else "still_evolving",
                    "n_snapshots": case.times_s.size,
                    "t_start_s": float(case.times_s[0]),
                    "t_final_s": float(case.times_s[-1]),
                    "previous_window_start_s": float(case.times_s[prev_slice.start]),
                    "final_window_start_s": float(case.times_s[final_slice.start]),
                    "relative_window_change_l2": window_change,
                    "relative_last_step_change_l2": last_change,
                    "final_window_integral": integral,
                    "final_window_peak": float(np.nanmax(final)),
                    "final_window_peak_z_m": float(case.z_m[peak_idx]),
                    "tolerance": tolerance,
                }
            )
    return rows


def plot_density_maps(case: SteadyCase, out_dir: Path, smooth_window: int) -> None:
    species = [spp for spp in case.species if "n_m" in case.mesh.get(spp, {})]
    if not species:
        return
    fig, axes = plt.subplots(len(species), 1, figsize=(13.0, 4.8 * len(species)), sharex=True, squeeze=False)
    for ax, spp in zip(axes.reshape(-1), species):
        values = smooth_z(case.mesh[spp]["n_m"], smooth_window)
        vmax = float(np.nanpercentile(values, 99.5)) if np.any(np.isfinite(values)) else 1.0
        im = ax.pcolormesh(case.z_m, case.times_s * 1.0e3, values.T, shading="auto", cmap="viridis", vmin=0.0, vmax=vmax)
        ax.set_ylabel("time [ms]")
        ax.set_title(f"{case.name}: {SPECIES_LABEL.get(spp, spp)} density")
        ax.grid(True, alpha=0.18)
        cbar = fig.colorbar(im, ax=ax, pad=0.01)
        cbar.set_label("n [m^-3]")
    axes[-1, 0].set_xlabel("z [m]")
    fig.tight_layout()
    fig.savefig(out_dir / f"{case.name}_density_time_z.png", dpi=220)
    plt.close(fig)


def plot_profile_windows(case: SteadyCase, out_dir: Path, smooth_window: int, window_fraction: float) -> None:
    if case.times_s.size < 2:
        return
    if case.times_s.size >= 4:
        prev_slice, final_slice = window_indices(case.times_s.size, window_fraction)
    else:
        prev_slice = slice(max(0, case.times_s.size - 2), max(1, case.times_s.size - 1))
        final_slice = slice(case.times_s.size - 1, case.times_s.size)

    variables = ("n_m", "Tpar_m", "Tper_m")
    fig, axes = plt.subplots(len(variables), len(case.species), figsize=(6.6 * len(case.species), 4.8 * len(variables)), squeeze=False)
    for col, spp in enumerate(case.species):
        for row, var in enumerate(variables):
            ax = axes[row, col]
            if var not in case.mesh.get(spp, {}):
                ax.axis("off")
                continue
            values = case.mesh[spp][var]
            prev = smooth_z(finite_mean_profile(values, prev_slice.start, prev_slice.stop), smooth_window)
            final = smooth_z(finite_mean_profile(values, final_slice.start, final_slice.stop), smooth_window)
            ax.plot(case.z_m, prev, color="0.45", lw=2.0, label="previous window")
            ax.plot(case.z_m, final, color="#d62728", lw=2.2, label="final window")
            ax.set_title(f"{SPECIES_LABEL.get(spp, spp)} {var}")
            ax.set_xlabel("z [m]")
            ax.set_ylabel("n [m^-3]" if var == "n_m" else "eV" if var.startswith("T") else "m/s")
            ax.grid(True, alpha=0.22)
            if row == 0 and col == 0:
                ax.legend(frameon=False)
    fig.suptitle(f"{case.name}: previous versus final profile windows")
    fig.tight_layout()
    fig.savefig(out_dir / f"{case.name}_profile_window_compare.png", dpi=220)
    plt.close(fig)


def plot_global_histories(case: SteadyCase, out_dir: Path) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(12.5, 12.0), sharex=True)
    for spp in case.species:
        label = SPECIES_LABEL.get(spp, spp)
        if "n_m" in case.mesh.get(spp, {}):
            density = case.mesh[spp]["n_m"]
            axes[0].plot(case.times_s * 1.0e3, np.trapz(density, case.z_m, axis=0), lw=2.0, label=label)
            axes[1].plot(case.times_s * 1.0e3, np.nanmax(density, axis=0), lw=2.0, label=label)
        for var, style in (("Tpar_m", "-"), ("Tper_m", "--")):
            if var in case.mesh.get(spp, {}):
                temp = case.mesh[spp][var]
                density = case.mesh[spp].get("n_m")
                if density is not None:
                    weights = np.maximum(density, 0.0)
                    mean_temp = np.sum(temp * weights, axis=0) / np.maximum(np.sum(weights, axis=0), 1.0e-300)
                else:
                    mean_temp = np.nanmean(temp, axis=0)
                axes[2].plot(case.times_s * 1.0e3, mean_temp, style, lw=2.0, label=f"{label} {var}")
    axes[0].set_ylabel("line integral n dz [m^-2]")
    axes[1].set_ylabel("peak n [m^-3]")
    axes[2].set_ylabel("density-weighted T [eV]")
    axes[2].set_xlabel("time [ms]")
    for ax in axes:
        ax.grid(True, alpha=0.22)
        ax.legend(frameon=False)
    fig.suptitle(f"{case.name}: global steady-state indicators")
    fig.tight_layout()
    fig.savefig(out_dir / f"{case.name}_global_time_histories.png", dpi=220)
    plt.close(fig)


def write_table(rows: Iterable[dict[str, object]], out_dir: Path) -> tuple[Path, Path]:
    rows = list(rows)
    csv_path = out_dir / "steady_state_metrics.csv"
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    md_path = out_dir / "steady_state_summary.md"
    lines = ["# PICOS++ MPEX Steady-State Profile Check", ""]
    if not rows:
        lines.append("No cases were analyzed.")
    else:
        by_case = sorted({str(row.get("case", "")) for row in rows})
        for case in by_case:
            case_rows = [row for row in rows if str(row.get("case", "")) == case]
            statuses = {str(row.get("status", "")) for row in case_rows}
            if "not_enough_snapshots" in statuses:
                verdict = "not enough snapshots"
            elif statuses == {"steady_like"}:
                verdict = "steady-like by configured tolerances"
            else:
                verdict = "still evolving by configured tolerances"
            lines.extend([f"## {case}", "", f"Verdict: **{verdict}**.", ""])
            lines.append("| species | variable | window change | last-step change | peak z [m] | final time [ms] | status |")
            lines.append("|---|---|---:|---:|---:|---:|---|")
            for row in case_rows:
                if row.get("status") == "not_enough_snapshots":
                    lines.append(f"|  | {row.get('quantity', '')} |  |  |  |  | {row.get('status')} |")
                    continue
                lines.append(
                    "| {species} | {variable} | {win:.4g} | {last:.4g} | {z:.4g} | {tf:.4g} | {status} |".format(
                        species=row.get("species", ""),
                        variable=row.get("variable", ""),
                        win=float(row.get("relative_window_change_l2", math.nan)),
                        last=float(row.get("relative_last_step_change_l2", math.nan)),
                        z=float(row.get("final_window_peak_z_m", math.nan)),
                        tf=1.0e3 * float(row.get("t_final_s", math.nan)),
                        status=row.get("status", ""),
                    )
                )
            lines.append("")
    md_path.write_text("\n".join(lines) + "\n")
    return csv_path, md_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path, help="HDF5 dir, output root, or output tar.gz to analyze.")
    parser.add_argument("--out-dir", type=Path, default=Path("validation/mpex_steady_state_profile_check"))
    parser.add_argument("--extract-dir", type=Path, default=None)
    parser.add_argument("--case-filter", default="steady", help="Only analyze HDF5 paths containing this substring; use empty string for all.")
    parser.add_argument("--window-fraction", type=float, default=0.25, help="Fraction of snapshots in each late-time comparison window.")
    parser.add_argument("--tolerance", type=float, default=0.05, help="Relative L2 threshold for steady-like status.")
    parser.add_argument("--smooth-window", type=int, default=5, help="Visualization-only moving-average window in z cells.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    hdf5_dirs: list[Path] = []
    for path in args.paths:
        found = find_hdf5_dirs(path, args.extract_dir)
        hdf5_dirs.extend(found)
    hdf5_dirs = sorted(set(hdf5_dirs))
    if args.case_filter:
        hdf5_dirs = [path for path in hdf5_dirs if args.case_filter in str(path)]

    if not hdf5_dirs:
        print("No PICOS HDF5 output directories found in the provided path(s).")
        print("The NERSC setup tarball contains input decks/scripts only; analyze the completed run output tarball or HDF5 directory.")
        return 2

    rows: list[dict[str, object]] = []
    for hdf5_dir in hdf5_dirs:
        case = read_case(hdf5_dir)
        plot_density_maps(case, args.out_dir, args.smooth_window)
        plot_profile_windows(case, args.out_dir, args.smooth_window, args.window_fraction)
        plot_global_histories(case, args.out_dir)
        rows.extend(metric_rows(case, args.window_fraction, args.tolerance))

    csv_path, md_path = write_table(rows, args.out_dir)
    print(md_path)
    print(csv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
