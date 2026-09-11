#!/usr/bin/env python3
"""Plot MPEX PICOS++ densities in the older ProtoLite preview style.

The reference ProtoLite script smooths mesh density with a moving average and
plots n_m as a surface over time and axial position. This script does the same
for the MPEX steady run plus the three ECH restart cases, and also writes 2D
time-z maps and final-window lineouts.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_ANALYSIS_ROOT = Path("/Users/78k/Desktop/MPEX_ECH_runs_analysis_20260911/MPEX_ECH_runs")
DEFAULT_OUT_DIR = Path("/Users/78k/Desktop/MPEX_ECH_runs_analysis_20260911/density_plots_protolite_style")

CASE_LAYOUT = [
    (
        "steady_2ms",
        "2 ms no-RF",
        DEFAULT_ANALYSIS_ROOT
        / "PICOS_NERSC_MPEX_steady_restart_2ms"
        / "picosFILES/outputFiles/mpex_scenario14_ex8_steady_2ms_p262144_coll_mpexprof_nersc_nonrel/HDF5",
    ),
    (
        "left_resonance_70ghz",
        "left resonance, 70 GHz",
        DEFAULT_ANALYSIS_ROOT
        / "PICOS_NERSC_MPEX_steady_then_ech_restart2ms"
        / "PICOS_NERSC_MPEX_triplet_restart2ms"
        / "left_resonance_70ghz/run_restart2ms/picosFILES/outputFiles"
        / "mpex_scenario14_ex8_left_resonance_70ghz_target8_p262144_100us_coll_mpexprof_restart2ms_nersc_nonrel/HDF5",
    ),
    (
        "right_resonance_70ghz",
        "right resonance, 70 GHz",
        DEFAULT_ANALYSIS_ROOT
        / "PICOS_NERSC_MPEX_steady_then_ech_restart2ms"
        / "PICOS_NERSC_MPEX_triplet_restart2ms"
        / "right_resonance_70ghz/run_restart2ms/picosFILES/outputFiles"
        / "mpex_scenario14_ex8_right_resonance_70ghz_target8_p262144_100us_coll_mpexprof_restart2ms_nersc_nonrel/HDF5",
    ),
    (
        "well_min_65p588ghz",
        "well minimum, 65.588 GHz",
        DEFAULT_ANALYSIS_ROOT
        / "PICOS_NERSC_MPEX_steady_then_ech_restart2ms"
        / "PICOS_NERSC_MPEX_triplet_restart2ms"
        / "well_min_65p588ghz/run_restart2ms/picosFILES/outputFiles"
        / "mpex_scenario14_ex8_well_min_65p588ghz_target8_p262144_100us_coll_mpexprof_restart2ms_nersc_nonrel/HDF5",
    ),
]

SPECIES_LABELS = {
    "spp_1": "D+",
    "spp_2": "electron",
}


@dataclass
class DensityCase:
    name: str
    label: str
    hdf5_dir: Path
    x_m: np.ndarray
    times_s: np.ndarray
    density: Dict[str, np.ndarray]
    b_m: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot PICOS++ MPEX density profiles.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--smooth-window", type=int, default=10)
    parser.add_argument("--final-window", type=int, default=10)
    return parser.parse_args()


def numeric_snapshot_names(h5: h5py.File) -> List[str]:
    return sorted([key for key in h5.keys() if key.isdigit()], key=lambda item: int(item))


def read_case(name: str, label: str, hdf5_dir: Path) -> DensityCase:
    if not hdf5_dir.exists():
        raise FileNotFoundError(hdf5_dir)
    with h5py.File(hdf5_dir / "main.h5", "r") as h5:
        x_m = np.asarray(h5["geometry/x_m"], dtype=np.float64).reshape(-1)

    p0 = hdf5_dir / "PARTICLES_FILE_0.h5"
    density: Dict[str, np.ndarray] = {}
    with h5py.File(p0, "r") as h5:
        snapshots = numeric_snapshot_names(h5)
        times = []
        for snapshot in snapshots:
            if f"{snapshot}/time" in h5:
                times.append(float(np.asarray(h5[f"{snapshot}/time"]).reshape(-1)[0]))
            else:
                times.append(float("nan"))
        for species in SPECIES_LABELS:
            values = []
            for snapshot in snapshots:
                dataset = f"{snapshot}/ions/{species}/n_m"
                if dataset not in h5:
                    raise KeyError(f"Missing {dataset} in {p0}")
                values.append(np.asarray(h5[dataset], dtype=np.float64).reshape(-1))
            density[species] = np.column_stack(values)

    with h5py.File(hdf5_dir / "FIELDS_FILE_0.h5", "r") as h5:
        b_m = np.asarray(h5[f"{snapshots[-1]}/fields/BX_m/x"], dtype=np.float64).reshape(-1)

    n = min(len(x_m), len(b_m), *(arr.shape[0] for arr in density.values()))
    return DensityCase(
        name=name,
        label=label,
        hdf5_dir=hdf5_dir,
        x_m=x_m[:n],
        times_s=np.asarray(times, dtype=np.float64),
        density={key: value[:n, :] for key, value in density.items()},
        b_m=b_m[:n],
    )


def moving_average_axis0(values: np.ndarray, window: int) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    if window <= 1:
        return values
    kernel = np.ones(window, dtype=np.float64) / float(window)
    return np.apply_along_axis(lambda column: np.convolve(column, kernel, mode="same"), 0, values)


def fixed_axes_3d(ax: plt.Axes, zmax: float) -> None:
    ax.set_zlim(0.0, zmax)
    ax.view_init(elev=24.0, azim=-135.0)
    ax.set_xlabel("time [us]", labelpad=8)
    ax.set_ylabel("z [m]", labelpad=8)
    ax.set_zlabel("n [m^-3]", labelpad=8)


def plot_density_surface(
    cases: Sequence[DensityCase],
    species: str,
    out_dir: Path,
    smooth_window: int,
) -> None:
    label = SPECIES_LABELS[species]
    smoothed = [moving_average_axis0(case.density[species], smooth_window) for case in cases]
    all_values = np.concatenate([arr[np.isfinite(arr)] for arr in smoothed if np.any(np.isfinite(arr))])
    zmax = float(np.nanpercentile(all_values, 99.7)) if all_values.size else 1.0
    zmax = max(zmax, 1.0)

    fig = plt.figure(figsize=(17.5, 12.5))
    for idx, (case, values) in enumerate(zip(cases, smoothed), start=1):
        ax = fig.add_subplot(2, 2, idx, projection="3d")
        t_grid, x_grid = np.meshgrid(case.times_s * 1.0e6, case.x_m)
        surf = ax.plot_surface(t_grid, x_grid, values, cmap="viridis", linewidth=0, antialiased=False)
        surf.set_clim(0.0, zmax)
        fixed_axes_3d(ax, zmax)
        ax.set_title(case.label)
    fig.suptitle(f"ProtoLite-style smoothed {label} density surfaces")
    fig.tight_layout()
    fig.savefig(out_dir / f"density_surface_{species}_protolite_style.png", dpi=220)
    plt.close(fig)


def plot_density_maps(
    cases: Sequence[DensityCase],
    species: str,
    out_dir: Path,
    smooth_window: int,
) -> None:
    label = SPECIES_LABELS[species]
    smoothed = [moving_average_axis0(case.density[species], smooth_window) for case in cases]
    all_values = np.concatenate([arr[np.isfinite(arr)] for arr in smoothed if np.any(np.isfinite(arr))])
    vmax = float(np.nanpercentile(all_values, 99.7)) if all_values.size else 1.0
    vmax = max(vmax, 1.0)

    fig, axes = plt.subplots(2, 2, figsize=(15.5, 10.0), sharex=True)
    axes = axes.reshape(-1)
    im = None
    for ax, case, values in zip(axes, cases, smoothed):
        im = ax.pcolormesh(case.x_m, case.times_s * 1.0e6, values.T, shading="auto", cmap="viridis", vmin=0.0, vmax=vmax)
        ax.set_title(case.label)
        ax.set_xlim(-2.0, 8.0)
        ax.set_ylabel("time [us]")
        ax.grid(True, alpha=0.18)
    for ax in axes[-2:]:
        ax.set_xlabel("z [m]")
    if im is not None:
        cbar = fig.colorbar(im, ax=list(axes), pad=0.02)
        cbar.set_label(f"{label} density [m^-3]")
    fig.suptitle(f"Smoothed {label} density time-z maps")
    fig.tight_layout(rect=(0.0, 0.0, 0.92, 0.95))
    fig.savefig(out_dir / f"density_map_{species}_protolite_style.png", dpi=220)
    plt.close(fig)


def plot_final_window_lineouts(
    cases: Sequence[DensityCase],
    out_dir: Path,
    smooth_window: int,
    final_window: int,
) -> List[Dict[str, object]]:
    fig, axes = plt.subplots(2, 2, figsize=(15.0, 9.6), sharex=True)
    rows: List[Dict[str, object]] = []
    for ax, case in zip(axes.reshape(-1), cases):
        start = max(0, case.times_s.size - final_window)
        for species, label in SPECIES_LABELS.items():
            final_density = np.nanmean(case.density[species][:, start:], axis=1)
            final_density = moving_average_axis0(final_density[:, None], smooth_window).reshape(-1)
            ax.plot(case.x_m, final_density, linewidth=2.2, label=label)
            rows.append(
                {
                    "case": case.name,
                    "species": label,
                    "final_window_outputs": case.times_s.size - start,
                    "final_time_us": float(case.times_s[-1] * 1.0e6),
                    "line_integrated_density": float(np.trapz(final_density, case.x_m)),
                    "peak_density_m^-3": float(np.nanmax(final_density)),
                    "peak_z_m": float(case.x_m[int(np.nanargmax(final_density))]),
                }
            )
        ax2 = ax.twinx()
        ax2.plot(case.x_m, case.b_m, color="0.35", linewidth=1.5, alpha=0.7, label="B")
        ax2.set_ylabel("B [T]", color="0.35")
        ax2.tick_params(axis="y", colors="0.35")
        ax.set_title(case.label)
        ax.set_ylabel("n [m^-3]")
        ax.set_xlim(-2.0, 8.0)
        ax.grid(True, alpha=0.25)
        ax.legend(frameon=False, loc="upper left")
    for ax in axes[-1, :]:
        ax.set_xlabel("z [m]")
    fig.suptitle("Final-window density profiles, ProtoLite smoothing")
    fig.tight_layout()
    fig.savefig(out_dir / "density_final_window_lineouts_protolite_style.png", dpi=220)
    plt.close(fig)
    return rows


def plot_normalized_final_overlay(
    cases: Sequence[DensityCase],
    out_dir: Path,
    smooth_window: int,
    final_window: int,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.8), sharex=True, sharey=True)
    colors = ["0.25", "#1f77b4", "#d62728", "#2ca02c"]
    for case, color in zip(cases, colors):
        start = max(0, case.times_s.size - final_window)
        for ax, species in zip(axes, ("spp_1", "spp_2")):
            density = np.nanmean(case.density[species][:, start:], axis=1)
            density = moving_average_axis0(density[:, None], smooth_window).reshape(-1)
            peak = float(np.nanmax(density)) if np.nanmax(density) > 0.0 else 1.0
            ax.plot(case.x_m, density / peak, linewidth=2.1, color=color, label=case.label)
            ax.set_title(f"{SPECIES_LABELS[species]} density")
            ax.set_xlabel("z [m]")
            ax.set_xlim(-2.0, 8.0)
            ax.grid(True, alpha=0.25)
    axes[0].set_ylabel("n / max(n)")
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle("Normalized final-window density comparison")
    fig.tight_layout()
    fig.savefig(out_dir / "density_final_window_normalized_overlay.png", dpi=220)
    plt.close(fig)


def write_rows(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    if not rows:
        return
    keys: List[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    cases = [read_case(name, label, path) for name, label, path in CASE_LAYOUT]
    for species in SPECIES_LABELS:
        plot_density_surface(cases, species, out_dir, args.smooth_window)
        plot_density_maps(cases, species, out_dir, args.smooth_window)
    rows = plot_final_window_lineouts(cases, out_dir, args.smooth_window, args.final_window)
    plot_normalized_final_overlay(cases, out_dir, args.smooth_window, args.final_window)
    write_rows(out_dir / "density_final_window_summary.csv", rows)

    print(f"Wrote density plots to {out_dir}")


if __name__ == "__main__":
    main()
