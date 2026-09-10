#!/usr/bin/env python3
"""Analyze MPEX ECH triplet PICOS++ outputs and plot electron EEDFs.

This script expects a PICOS_NERSC_MPEX_triplet-style run directory with cases
such as left_resonance_70ghz, right_resonance_70ghz, and well_min_65p588ghz.
It aggregates all PARTICLES_FILE_*.h5 files for the kinetic electron species
and writes weighted EEDF and energy-versus-z diagnostics.
"""

from __future__ import annotations

import argparse
import csv
import math
import tarfile
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


E_CHARGE = 1.602176634e-19
M_E = 9.1093837015e-31

CASE_LABELS = {
    "left_resonance_70ghz": "left resonance, 70 GHz",
    "right_resonance_70ghz": "right resonance, 70 GHz",
    "well_min_65p588ghz": "well minimum, 65.588 GHz",
}

CASE_COLORS = {
    "left_resonance_70ghz": "#1f77b4",
    "right_resonance_70ghz": "#d62728",
    "well_min_65p588ghz": "#2ca02c",
}

SELECTED_TIME_FRACTIONS = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)


@dataclass
class CaseInfo:
    name: str
    label: str
    hdf5_dir: Path
    input_file: Optional[Path]
    tag: str
    particles: List[Path]
    fields_file: Path
    main_file: Path
    input_values: Dict[str, str]
    snapshots: List[str]
    times_s: List[float]
    resonance_z: Optional[float]
    rf_frequency_hz: Optional[float]


@dataclass
class PhaseSpaceData:
    z: np.ndarray
    energy_ev: np.ndarray
    weights: np.ndarray
    v_parallel: np.ndarray
    v_perp: np.ndarray
    b_local: np.ndarray
    phi_local: np.ndarray
    mu_j_per_t: np.ndarray
    pitch_sin2: np.ndarray
    loss_cone_sin2: np.ndarray
    trapped_magnetic: np.ndarray
    trapped_effective: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze PICOS++ MPEX ECH triplet HDF5 outputs and plot EEDFs."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=None,
        help="Extracted PICOS_NERSC_MPEX_triplet_mpexprof directory.",
    )
    parser.add_argument(
        "--archive",
        type=Path,
        default=None,
        help="Optional tar.gz archive to extract before analysis.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("PICOS_NERSC_MPEX_triplet_mpexprof_analysis"),
        help="Directory where plots and CSV summaries are written.",
    )
    parser.add_argument(
        "--electron-species",
        default="spp_2",
        help="Kinetic electron species group name under /ions.",
    )
    parser.add_argument(
        "--energy-max",
        type=float,
        default=None,
        help="Maximum energy in eV for EEDF and E-z plots. Default: inferred.",
    )
    parser.add_argument(
        "--energy-bins",
        type=int,
        default=300,
        help="Number of linear energy bins for EEDF plots.",
    )
    parser.add_argument(
        "--z-bins",
        type=int,
        default=260,
        help="Number of z bins for E-z histograms.",
    )
    parser.add_argument(
        "--fig17-energy-max",
        type=float,
        default=3500.0,
        help="Maximum energy in eV for individual Fig. 17-style plots.",
    )
    parser.add_argument(
        "--figure-current-a",
        type=float,
        default=1900.0,
        help="Current label in amperes for individual Fig. 17-style plot titles.",
    )
    parser.add_argument(
        "--hot-thresholds",
        default="100,300,500",
        help="Comma-separated electron energy thresholds in eV for hot-population plots.",
    )
    return parser.parse_args()


def extract_archive(archive: Path) -> Path:
    if not archive.exists():
        raise FileNotFoundError(f"Archive does not exist: {archive}")
    tmp_root = Path(tempfile.mkdtemp(prefix="picos_mpexprof_eedf_"))
    with tarfile.open(archive, "r:gz") as tf:
        tf.extractall(tmp_root)
    roots = [p for p in tmp_root.iterdir() if p.is_dir()]
    if len(roots) == 1:
        return roots[0]
    return tmp_root


def numeric_suffix(path: Path) -> int:
    stem = path.stem
    try:
        return int(stem.rsplit("_", 1)[1])
    except (IndexError, ValueError):
        return 10**9


def parse_input_file(path: Optional[Path]) -> Dict[str, str]:
    values: Dict[str, str] = {}
    if path is None or not path.exists():
        return values
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.split("//", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 2:
            values[parts[0]] = parts[1]
    return values


def as_float(values: Dict[str, str], key: str) -> Optional[float]:
    try:
        return float(values[key])
    except (KeyError, ValueError):
        return None


def find_input_file(picos_files_dir: Path, tag: str) -> Optional[Path]:
    input_dir = picos_files_dir / "inputFiles"
    if not input_dir.exists():
        return None
    exact = input_dir / f"input_file_{tag}.input"
    if exact.exists():
        return exact
    matches = sorted(input_dir.glob("input_file_*.input"))
    return matches[0] if matches else None


def hdf5_snapshots_and_times(particle_file: Path) -> Tuple[List[str], List[float]]:
    with h5py.File(particle_file, "r") as h5:
        snapshots = sorted([k for k in h5.keys() if k.isdigit()], key=lambda x: int(x))
        times: List[float] = []
        for snap in snapshots:
            if "time" in h5[snap]:
                times.append(float(np.asarray(h5[snap]["time"]).reshape(-1)[0]))
            elif "time[s]" in h5[snap].attrs:
                times.append(float(np.asarray(h5[snap].attrs["time[s]"]).reshape(-1)[0]))
            else:
                times.append(float("nan"))
    return snapshots, times


def discover_cases(root: Path) -> List[CaseInfo]:
    hdf5_dirs = sorted(root.glob("*/run_*interactive/picosFILES/outputFiles/*/HDF5"))
    if not hdf5_dirs:
        hdf5_dirs = sorted(root.glob("*/run*/picosFILES/outputFiles/*/HDF5"))
    if not hdf5_dirs:
        raise FileNotFoundError(f"No PICOS HDF5 output directories found under {root}")

    cases: List[CaseInfo] = []
    for hdf5_dir in hdf5_dirs:
        tag = hdf5_dir.parent.name
        picos_files_dir = hdf5_dir.parents[2]
        case_dir = hdf5_dir.parents[4]
        case_name = case_dir.name
        particle_files = sorted(hdf5_dir.glob("PARTICLES_FILE_*.h5"), key=numeric_suffix)
        fields_file = hdf5_dir / "FIELDS_FILE_0.h5"
        main_file = hdf5_dir / "main.h5"
        if not particle_files:
            raise FileNotFoundError(f"No PARTICLES_FILE_*.h5 in {hdf5_dir}")
        if not fields_file.exists():
            raise FileNotFoundError(f"Missing fields file: {fields_file}")
        if not main_file.exists():
            raise FileNotFoundError(f"Missing main file: {main_file}")
        input_file = find_input_file(picos_files_dir, tag)
        input_values = parse_input_file(input_file)
        snapshots, times_s = hdf5_snapshots_and_times(particle_files[0])
        x1 = as_float(input_values, "RF_electron_x1")
        x2 = as_float(input_values, "RF_electron_x2")
        resonance_z = None
        if x1 is not None and x2 is not None:
            # These NERSC cases were generated with x1 = z_res - 0.02 m
            # and x2 = target position. Store the resonance estimate for plots.
            resonance_z = x1 + 0.02
        cases.append(
            CaseInfo(
                name=case_name,
                label=CASE_LABELS.get(case_name, case_name),
                hdf5_dir=hdf5_dir,
                input_file=input_file,
                tag=tag,
                particles=particle_files,
                fields_file=fields_file,
                main_file=main_file,
                input_values=input_values,
                snapshots=snapshots,
                times_s=times_s,
                resonance_z=resonance_z,
                rf_frequency_hz=as_float(input_values, "RF_electron_freq"),
            )
        )
    return sorted(cases, key=lambda c: c.name)


def read_particles(
    case: CaseInfo, snapshot: str, electron_species: str
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    z_chunks: List[np.ndarray] = []
    e_chunks: List[np.ndarray] = []
    w_chunks: List[np.ndarray] = []
    base = f"{snapshot}/ions/{electron_species}"
    for path in case.particles:
        with h5py.File(path, "r") as h5:
            if base not in h5:
                raise KeyError(f"Missing group {base} in {path}")
            grp = h5[base]
            x = np.asarray(grp["X_p"], dtype=np.float64).reshape(-1)
            weight = np.asarray(grp["a_p"], dtype=np.float64).reshape(-1)
            velocity = np.asarray(grp["V_p"], dtype=np.float64)
            if velocity.shape[0] != 2:
                raise ValueError(f"Expected V_p shape [2, n], got {velocity.shape} in {path}")
            energy_ev = 0.5 * M_E * (velocity[0] ** 2 + velocity[1] ** 2) / E_CHARGE
            valid = np.isfinite(x) & np.isfinite(energy_ev) & np.isfinite(weight) & (weight > 0)
            z_chunks.append(x[valid])
            e_chunks.append(energy_ev[valid])
            w_chunks.append(weight[valid])
    return np.concatenate(z_chunks), np.concatenate(e_chunks), np.concatenate(w_chunks)


def effective_barriers_for_particles(
    b_grid: np.ndarray,
    phi_grid: np.ndarray,
    particle_indices: np.ndarray,
    mu_j_per_t: np.ndarray,
    charge_c: float,
    chunk_size: int = 4096,
) -> Tuple[np.ndarray, np.ndarray]:
    left_barrier = np.empty(mu_j_per_t.size, dtype=np.float64)
    right_barrier = np.empty(mu_j_per_t.size, dtype=np.float64)
    grid_index = np.arange(b_grid.size, dtype=np.int32)[None, :]
    for start in range(0, mu_j_per_t.size, chunk_size):
        stop = min(start + chunk_size, mu_j_per_t.size)
        mu_chunk = mu_j_per_t[start:stop]
        idx_chunk = particle_indices[start:stop][:, None]
        effective_potential = mu_chunk[:, None] * b_grid[None, :] + charge_c * phi_grid[None, :]
        left_barrier[start:stop] = np.max(
            np.where(grid_index <= idx_chunk, effective_potential, -np.inf), axis=1
        )
        right_barrier[start:stop] = np.max(
            np.where(grid_index >= idx_chunk, effective_potential, -np.inf), axis=1
        )
    return left_barrier, right_barrier


def read_phase_space(
    case: CaseInfo, snapshot: str, electron_species: str
) -> PhaseSpaceData:
    z_chunks: List[np.ndarray] = []
    e_chunks: List[np.ndarray] = []
    w_chunks: List[np.ndarray] = []
    vpar_chunks: List[np.ndarray] = []
    vperp_chunks: List[np.ndarray] = []
    base = f"{snapshot}/ions/{electron_species}"
    for path in case.particles:
        with h5py.File(path, "r") as h5:
            if base not in h5:
                raise KeyError(f"Missing group {base} in {path}")
            grp = h5[base]
            x = np.asarray(grp["X_p"], dtype=np.float64).reshape(-1)
            weight = np.asarray(grp["a_p"], dtype=np.float64).reshape(-1)
            velocity = np.asarray(grp["V_p"], dtype=np.float64)
            if velocity.shape[0] != 2:
                raise ValueError(f"Expected V_p shape [2, n], got {velocity.shape} in {path}")
            v_parallel = velocity[0]
            v_perp = velocity[1]
            v2 = v_parallel**2 + v_perp**2
            energy_ev = 0.5 * M_E * v2 / E_CHARGE
            valid = (
                np.isfinite(x)
                & np.isfinite(energy_ev)
                & np.isfinite(weight)
                & np.isfinite(v_parallel)
                & np.isfinite(v_perp)
                & (weight > 0)
            )
            z_chunks.append(x[valid])
            e_chunks.append(energy_ev[valid])
            w_chunks.append(weight[valid])
            vpar_chunks.append(v_parallel[valid])
            vperp_chunks.append(v_perp[valid])

    z = np.concatenate(z_chunks)
    energy_ev = np.concatenate(e_chunks)
    weights = np.concatenate(w_chunks)
    v_parallel = np.concatenate(vpar_chunks)
    v_perp = np.concatenate(vperp_chunks)

    x_m, b_m, _, phi_m = read_field_profile(case, snapshot)
    order = np.argsort(x_m)
    x_grid = x_m[order]
    b_grid = b_m[order]
    phi_grid = phi_m[order]
    b_local = np.interp(z, x_grid, b_grid, left=b_grid[0], right=b_grid[-1])
    phi_local = np.interp(z, x_grid, phi_grid, left=phi_grid[0], right=phi_grid[-1])

    v2 = v_parallel**2 + v_perp**2
    vperp2 = v_perp**2
    with np.errstate(divide="ignore", invalid="ignore"):
        pitch_sin2 = np.where(v2 > 0.0, vperp2 / v2, 0.0)
        mu_j_per_t = np.where(b_local > 0.0, 0.5 * M_E * vperp2 / b_local, 0.0)

    b_left_max_grid = np.maximum.accumulate(b_grid)
    b_right_max_grid = np.maximum.accumulate(b_grid[::-1])[::-1]
    b_left_max = np.interp(z, x_grid, b_left_max_grid, left=b_left_max_grid[0], right=b_left_max_grid[-1])
    b_right_max = np.interp(z, x_grid, b_right_max_grid, left=b_right_max_grid[0], right=b_right_max_grid[-1])
    with np.errstate(divide="ignore", invalid="ignore"):
        loss_cone_sin2 = np.maximum(b_local / b_left_max, b_local / b_right_max)
    loss_cone_sin2 = np.clip(loss_cone_sin2, 0.0, 1.0)
    trapped_magnetic = pitch_sin2 >= loss_cone_sin2

    particle_indices = np.clip(np.searchsorted(x_grid, z), 0, x_grid.size - 1)
    total_energy_j = energy_ev * E_CHARGE
    electron_hamiltonian = total_energy_j - E_CHARGE * phi_local
    left_u_max, right_u_max = effective_barriers_for_particles(
        b_grid, phi_grid, particle_indices, mu_j_per_t, -E_CHARGE
    )
    trapped_effective = (left_u_max > electron_hamiltonian) & (right_u_max > electron_hamiltonian)

    return PhaseSpaceData(
        z=z,
        energy_ev=energy_ev,
        weights=weights,
        v_parallel=v_parallel,
        v_perp=v_perp,
        b_local=b_local,
        phi_local=phi_local,
        mu_j_per_t=mu_j_per_t,
        pitch_sin2=pitch_sin2,
        loss_cone_sin2=loss_cone_sin2,
        trapped_magnetic=trapped_magnetic,
        trapped_effective=trapped_effective,
    )


def read_field_profile(case: CaseInfo, snapshot: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(case.main_file, "r") as h5:
        x_m = np.asarray(h5["geometry/x_m"], dtype=np.float64).reshape(-1)
    with h5py.File(case.fields_file, "r") as h5:
        grp = h5[snapshot]
        b = np.asarray(grp["fields/BX_m/x"], dtype=np.float64).reshape(-1)
        e = np.asarray(grp["fields/EX_m/x"], dtype=np.float64).reshape(-1)
        phi = np.asarray(grp["fields/Phi_m/x"], dtype=np.float64).reshape(-1)
    n = min(len(x_m), len(b), len(e), len(phi))
    return x_m[:n], b[:n], e[:n], phi[:n]


def weighted_percentile(values: np.ndarray, weights: np.ndarray, percentile: float) -> float:
    if values.size == 0:
        return float("nan")
    sorter = np.argsort(values)
    vals = values[sorter]
    w = weights[sorter]
    cdf = np.cumsum(w)
    if cdf[-1] <= 0:
        return float(np.percentile(values, percentile))
    target = percentile / 100.0 * cdf[-1]
    return float(vals[np.searchsorted(cdf, target, side="left")])


def weighted_fraction_above(values: np.ndarray, weights: np.ndarray, threshold: float) -> float:
    total = float(np.sum(weights))
    if total <= 0:
        return float("nan")
    return float(np.sum(weights[values >= threshold]) / total)


def summarize_snapshot(
    case: CaseInfo, snap_index: int, z: np.ndarray, energy: np.ndarray, weights: np.ndarray
) -> Dict[str, float | str | int]:
    total_weight = float(np.sum(weights))
    mean_e = float(np.average(energy, weights=weights)) if total_weight > 0 else float("nan")
    z_mean = float(np.average(z, weights=weights)) if total_weight > 0 else float("nan")
    return {
        "case": case.name,
        "label": case.label,
        "snapshot": int(case.snapshots[snap_index]),
        "time_s": case.times_s[snap_index],
        "particle_count": int(energy.size),
        "total_weight": total_weight,
        "mean_energy_eV": mean_e,
        "median_energy_eV": weighted_percentile(energy, weights, 50.0),
        "p90_energy_eV": weighted_percentile(energy, weights, 90.0),
        "p99_energy_eV": weighted_percentile(energy, weights, 99.0),
        "p999_energy_eV": weighted_percentile(energy, weights, 99.9),
        "max_energy_eV": float(np.max(energy)) if energy.size else float("nan"),
        "frac_weight_E_gt_100eV": weighted_fraction_above(energy, weights, 100.0),
        "frac_weight_E_gt_300eV": weighted_fraction_above(energy, weights, 300.0),
        "frac_weight_E_gt_500eV": weighted_fraction_above(energy, weights, 500.0),
        "frac_weight_E_gt_1000eV": weighted_fraction_above(energy, weights, 1000.0),
        "weighted_mean_z_m": z_mean,
        "weighted_z_p10_m": weighted_percentile(z, weights, 10.0),
        "weighted_z_p90_m": weighted_percentile(z, weights, 90.0),
        "resonance_z_m": case.resonance_z if case.resonance_z is not None else float("nan"),
        "rf_frequency_hz": case.rf_frequency_hz if case.rf_frequency_hz is not None else float("nan"),
    }


def write_rows(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        return
    keys: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def selected_snapshot_indices(count: int) -> List[int]:
    indices = sorted(set(int(round(frac * (count - 1))) for frac in SELECTED_TIME_FRACTIONS))
    return indices


def eedf_histogram(
    energy: np.ndarray, weights: np.ndarray, bins: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    hist, edges = np.histogram(energy, bins=bins, weights=weights)
    width = np.diff(edges)
    total = np.sum(hist)
    pdf = hist / (total * width) if total > 0 else np.zeros_like(hist)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, pdf


def set_axis_style(ax: plt.Axes) -> None:
    ax.grid(True, alpha=0.25, linewidth=0.8)
    ax.tick_params(direction="out")


def plot_final_eedf(
    cases: Sequence[CaseInfo],
    final_data: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]],
    bins: np.ndarray,
    electron_species: str,
    out_dir: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 6.5))
    for case in cases:
        _, energy, weights = final_data[case.name]
        centers, pdf = eedf_histogram(energy, weights, bins)
        ax.semilogy(
            centers,
            np.maximum(pdf, 1e-30),
            label=case.label,
            color=CASE_COLORS.get(case.name),
            linewidth=2.2,
        )
    z0, e0, w0 = read_particles(cases[0], cases[0].snapshots[0], electron_species)
    centers0, pdf0 = eedf_histogram(e0, w0, bins)
    ax.semilogy(centers0, np.maximum(pdf0, 1e-30), "--", color="0.35", label="initial")
    ax.set_xlabel("electron kinetic energy [eV]")
    ax.set_ylabel("weighted probability density [1/eV]")
    ax.set_title("Global electron energy distribution at final output")
    ax.set_xlim(bins[0], bins[-1])
    ax.legend(frameon=False)
    set_axis_style(ax)
    fig.tight_layout()
    fig.savefig(out_dir / "eedf_final_global.png", dpi=220)
    plt.close(fig)


def plot_final_eedf_logbins(
    cases: Sequence[CaseInfo],
    final_data: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]],
    energy_max: float,
    electron_species: str,
    out_dir: Path,
) -> None:
    bins = np.geomspace(0.2, max(energy_max, 1.0), 180)
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 6.0))
    for case in cases:
        _, energy, weights = final_data[case.name]
        centers, pdf = eedf_histogram(energy, weights, bins)
        axes[0].loglog(
            centers,
            np.maximum(pdf, 1e-30),
            label=case.label,
            color=CASE_COLORS.get(case.name),
            linewidth=2.2,
        )

        thresholds = np.geomspace(1.0, max(energy_max, 1.0), 220)
        total_weight = np.sum(weights)
        tail = np.array(
            [
                np.sum(weights[energy >= threshold]) / total_weight if total_weight > 0 else np.nan
                for threshold in thresholds
            ]
        )
        axes[1].loglog(
            thresholds,
            np.maximum(tail, 1e-12),
            label=case.label,
            color=CASE_COLORS.get(case.name),
            linewidth=2.2,
        )

    _, e0, w0 = read_particles(cases[0], cases[0].snapshots[0], electron_species)
    centers0, pdf0 = eedf_histogram(e0, w0, bins)
    axes[0].loglog(
        centers0, np.maximum(pdf0, 1e-30), "--", color="0.35", linewidth=2.0, label="initial"
    )

    axes[0].set_xlabel("electron kinetic energy [eV]")
    axes[0].set_ylabel("weighted probability density [1/eV]")
    axes[0].set_title("Log-binned EEDF")
    axes[1].set_xlabel("threshold energy [eV]")
    axes[1].set_ylabel("weighted fraction above threshold")
    axes[1].set_title("Cumulative high-energy tail")
    for ax in axes:
        ax.set_xlim(0.2, energy_max)
        ax.legend(frameon=False, fontsize=9)
        set_axis_style(ax)
    fig.suptitle("Final electron energy distribution with tail diagnostics")
    fig.tight_layout()
    fig.savefig(out_dir / "eedf_final_logbins_and_tail.png", dpi=220)
    plt.close(fig)


def plot_time_evolution_eedf(
    cases: Sequence[CaseInfo],
    electron_species: str,
    bins: np.ndarray,
    out_dir: Path,
) -> None:
    for case in cases:
        fig, ax = plt.subplots(figsize=(10.5, 6.5))
        for idx in selected_snapshot_indices(len(case.snapshots)):
            z, energy, weights = read_particles(case, case.snapshots[idx], electron_species)
            centers, pdf = eedf_histogram(energy, weights, bins)
            label = f"{case.times_s[idx] * 1e6:.1f} us"
            ax.semilogy(centers, np.maximum(pdf, 1e-30), linewidth=1.8, label=label)
        ax.set_xlabel("electron kinetic energy [eV]")
        ax.set_ylabel("weighted probability density [1/eV]")
        ax.set_title(f"EEDF time evolution: {case.label}")
        ax.set_xlim(bins[0], bins[-1])
        ax.legend(frameon=False, ncol=2)
        set_axis_style(ax)
        fig.tight_layout()
        fig.savefig(out_dir / f"eedf_time_evolution_{case.name}.png", dpi=220)
        plt.close(fig)


def plot_region_eedf(
    cases: Sequence[CaseInfo],
    final_data: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]],
    bins: np.ndarray,
    out_dir: Path,
) -> None:
    fig, axes = plt.subplots(len(cases), 1, figsize=(10.5, 3.7 * len(cases)), sharex=True)
    if len(cases) == 1:
        axes = [axes]
    for ax, case in zip(axes, cases):
        z, energy, weights = final_data[case.name]
        z_res = case.resonance_z if case.resonance_z is not None else 3.0
        regions = [
            ("all", np.ones_like(z, dtype=bool)),
            ("upstream", z < z_res - 0.5),
            ("near RF/resonance", np.abs(z - z_res) <= 0.25),
            ("downstream", (z > z_res + 0.5) & (z < 7.0)),
            ("target side", z >= 7.0),
        ]
        for label, mask in regions:
            if np.count_nonzero(mask) < 10:
                continue
            centers, pdf = eedf_histogram(energy[mask], weights[mask], bins)
            ax.semilogy(centers, np.maximum(pdf, 1e-30), linewidth=1.8, label=label)
        ax.set_ylabel("PDF [1/eV]")
        ax.set_title(case.label)
        ax.legend(frameon=False, fontsize=9, ncol=3)
        set_axis_style(ax)
    axes[-1].set_xlabel("electron kinetic energy [eV]")
    axes[-1].set_xlim(bins[0], bins[-1])
    fig.suptitle("Final EEDF by axial region", y=0.995)
    fig.tight_layout()
    fig.savefig(out_dir / "eedf_final_regions.png", dpi=220)
    plt.close(fig)


def plot_energy_metrics(time_rows: List[Dict[str, object]], out_dir: Path) -> None:
    metrics = [
        ("mean_energy_eV", "weighted mean [eV]"),
        ("p99_energy_eV", "weighted p99 [eV]"),
        ("p999_energy_eV", "weighted p99.9 [eV]"),
        ("max_energy_eV", "max [eV]"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.2), sharex=True)
    by_case: Dict[str, List[Dict[str, object]]] = {}
    for row in time_rows:
        by_case.setdefault(str(row["case"]), []).append(row)
    for ax, (metric, ylabel) in zip(axes.reshape(-1), metrics):
        for case_name, rows in sorted(by_case.items()):
            rows = sorted(rows, key=lambda r: float(r["time_s"]))
            t_us = np.array([float(r["time_s"]) * 1e6 for r in rows])
            y = np.array([float(r[metric]) for r in rows])
            ax.plot(
                t_us,
                y,
                label=CASE_LABELS.get(case_name, case_name),
                color=CASE_COLORS.get(case_name),
                linewidth=2.0,
            )
        ax.set_ylabel(ylabel)
        set_axis_style(ax)
    for ax in axes[-1, :]:
        ax.set_xlabel("time [us]")
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle("Electron energy metrics from weighted particle ensemble")
    fig.tight_layout()
    fig.savefig(out_dir / "energy_metrics_vs_time.png", dpi=220)
    plt.close(fig)


def plot_tail_fraction(time_rows: List[Dict[str, object]], out_dir: Path) -> None:
    thresholds = [
        ("frac_weight_E_gt_100eV", ">100 eV"),
        ("frac_weight_E_gt_300eV", ">300 eV"),
        ("frac_weight_E_gt_500eV", ">500 eV"),
        ("frac_weight_E_gt_1000eV", ">1 keV"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.2), sharex=True)
    by_case: Dict[str, List[Dict[str, object]]] = {}
    for row in time_rows:
        by_case.setdefault(str(row["case"]), []).append(row)
    for ax, (metric, title) in zip(axes.reshape(-1), thresholds):
        for case_name, rows in sorted(by_case.items()):
            rows = sorted(rows, key=lambda r: float(r["time_s"]))
            t_us = np.array([float(r["time_s"]) * 1e6 for r in rows])
            y = np.array([float(r[metric]) for r in rows])
            ax.semilogy(
                t_us,
                np.maximum(y, 1e-12),
                label=CASE_LABELS.get(case_name, case_name),
                color=CASE_COLORS.get(case_name),
                linewidth=2.0,
            )
        ax.set_title(title)
        ax.set_ylabel("weighted fraction")
        set_axis_style(ax)
    for ax in axes[-1, :]:
        ax.set_xlabel("time [us]")
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle("High-energy electron tail growth")
    fig.tight_layout()
    fig.savefig(out_dir / "tail_fraction_vs_time.png", dpi=220)
    plt.close(fig)


def plot_energy_vs_z(
    cases: Sequence[CaseInfo],
    final_data: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]],
    energy_max: float,
    z_bins: int,
    out_dir: Path,
) -> None:
    histograms = []
    extent_by_case = []
    for case in cases:
        z, energy, weights = final_data[case.name]
        z_edges = np.linspace(-2.0, 8.0, z_bins + 1)
        e_edges = np.linspace(0.0, energy_max, 260)
        hist, _, _ = np.histogram2d(z, energy, bins=(z_edges, e_edges), weights=weights)
        dz = np.diff(z_edges)[:, None]
        de = np.diff(e_edges)[None, :]
        total = np.sum(hist)
        density = hist / (total * dz * de) if total > 0 else hist
        log_density = np.full_like(density, np.nan, dtype=np.float64)
        positive = density > 0
        log_density[positive] = np.log10(density[positive])
        histograms.append(log_density)
        extent_by_case.append((z_edges[0], z_edges[-1], e_edges[0], e_edges[-1]))

    finite_values = np.concatenate([h[np.isfinite(h)] for h in histograms if np.any(np.isfinite(h))])
    vmin, vmax = np.percentile(finite_values, [5, 99.5]) if finite_values.size else (-8.0, -3.0)

    fig, axes = plt.subplots(1, len(cases), figsize=(7.2 * len(cases), 6.2), sharey=True)
    if len(cases) == 1:
        axes = [axes]
    im = None
    for ax, case, log_density, extent in zip(axes, cases, histograms, extent_by_case):
        im = ax.imshow(
            log_density.T,
            origin="lower",
            extent=extent,
            aspect="auto",
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
        )
        x_m, b_m, _, _ = read_field_profile(case, case.snapshots[-1])
        ax2 = ax.twinx()
        ax2.plot(x_m, b_m, color="#d62728", linewidth=2.0, alpha=0.9)
        ax2.set_ylim(0.0, max(1.1 * float(np.nanmax(b_m)), 1.5))
        ax2.set_ylabel("B [T]", color="#d62728")
        ax2.tick_params(axis="y", colors="#d62728")
        if case.resonance_z is not None:
            ax.axvline(case.resonance_z, color="white", linestyle="--", linewidth=1.8)
        ax.axvline(8.0, color="lime", linewidth=2.0)
        ax.set_xlabel("z [m]")
        ax.set_title(case.label)
        set_axis_style(ax)
    axes[0].set_ylabel("electron kinetic energy [eV]")
    if im is not None:
        cbar = fig.colorbar(im, ax=list(axes), pad=0.02)
        cbar.set_label("log10 weighted PDF [1/(m eV)]")
    fig.suptitle("Final electron energy distribution versus z")
    fig.tight_layout(rect=(0.0, 0.0, 0.93, 0.95))
    fig.savefig(out_dir / "electron_energy_vs_z_final.png", dpi=220)
    plt.close(fig)


def plot_fig17_style_maps(
    cases: Sequence[CaseInfo],
    final_data: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]],
    energy_max: float,
    current_a: float,
    out_dir: Path,
) -> None:
    rc = {
        "font.size": 24,
        "axes.titlesize": 28,
        "axes.labelsize": 26,
        "xtick.labelsize": 22,
        "ytick.labelsize": 22,
        "legend.fontsize": 18,
        "figure.titlesize": 30,
    }
    with plt.rc_context(rc):
        for case in cases:
            z, energy, weights = final_data[case.name]
            z_edges = np.linspace(-2.0, 8.0, 240)
            e_edges = np.linspace(0.0, energy_max, 220)
            hist, _, _ = np.histogram2d(z, energy, bins=(z_edges, e_edges), weights=weights)
            total = float(np.sum(hist))
            probability = hist / total if total > 0 else hist
            log_probability = np.full_like(probability, np.nan, dtype=np.float64)
            mask = probability > 0
            log_probability[mask] = np.log10(probability[mask])

            fig, ax = plt.subplots(figsize=(16.0, 8.2))
            im = ax.imshow(
                np.ma.masked_invalid(log_probability.T),
                origin="lower",
                extent=(z_edges[0], z_edges[-1], e_edges[0], e_edges[-1]),
                aspect="auto",
                cmap="viridis",
                vmin=-6.0,
                vmax=-2.0,
                interpolation="nearest",
            )

            if case.resonance_z is not None:
                ax.axvline(case.resonance_z, color="red", linestyle=":", linewidth=5.0)
            ax.axvline(8.0, color="limegreen", linewidth=5.0)

            x_m, b_m, _, _ = read_field_profile(case, case.snapshots[-1])
            ax2 = ax.twinx()
            ax2.plot(x_m, b_m, color="red", linewidth=4.0)
            ax2.set_ylabel("B0 [T]", color="red")
            ax2.tick_params(axis="y", colors="red", width=2.0, length=7)
            ax2.set_ylim(0.0, max(1.45, 1.05 * float(np.nanmax(b_m))))
            ax2.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

            ax.set_xlim(-2.0, 8.0)
            ax.set_ylim(0.0, energy_max)
            ax.set_xlabel("z [m]")
            ax.set_ylabel("E_e [eV]")
            ax.set_title(f"I: {current_a:.0f} [A] (PICOS++ nonrel)\n{case.label}")
            ax.grid(True, color="0.7", alpha=0.3, linewidth=1.4)
            ax.tick_params(width=2.0, length=7)

            cax = inset_axes(
                ax,
                width="4.0%",
                height="76%",
                loc="upper right",
                bbox_to_anchor=(-0.13, -0.02, 1.0, 1.0),
                bbox_transform=ax.transAxes,
                borderpad=0,
            )
            cbar = fig.colorbar(im, cax=cax)
            cbar.ax.set_title("log10(f(E,z))", fontsize=22, pad=12)
            cbar.set_ticks([-6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2])
            cbar.ax.tick_params(width=2.0, length=7)

            fig.subplots_adjust(left=0.10, right=0.86, bottom=0.14, top=0.82)
            fig.savefig(out_dir / f"fig17_style_{case.name}.png", dpi=220)
            plt.close(fig)


def plot_single_fig17_hot_map(
    case: CaseInfo,
    phase: PhaseSpaceData,
    hot_mask: np.ndarray,
    energy_max: float,
    current_a: float,
    threshold_ev: float,
    out_path: Path,
) -> None:
    rc = {
        "font.size": 24,
        "axes.titlesize": 22,
        "axes.labelsize": 26,
        "xtick.labelsize": 22,
        "ytick.labelsize": 22,
        "legend.fontsize": 18,
        "figure.titlesize": 30,
    }
    with plt.rc_context(rc):
        fig, ax = plt.subplots(figsize=(16.0, 8.2))
        z_edges = np.linspace(-2.0, 8.0, 240)
        e_edges = np.linspace(0.0, energy_max, 220)
        if np.count_nonzero(hot_mask) > 0 and np.sum(phase.weights[hot_mask]) > 0.0:
            hist, _, _ = np.histogram2d(
                phase.z[hot_mask],
                phase.energy_ev[hot_mask],
                bins=(z_edges, e_edges),
                weights=phase.weights[hot_mask],
            )
            probability = hist / np.sum(hist)
            log_probability = np.full_like(probability, np.nan, dtype=np.float64)
            positive = probability > 0.0
            log_probability[positive] = np.log10(probability[positive])
            image_data = np.ma.masked_invalid(log_probability.T)
        else:
            image_data = np.ma.masked_all((len(e_edges) - 1, len(z_edges) - 1))

        im = ax.imshow(
            image_data,
            origin="lower",
            extent=(z_edges[0], z_edges[-1], e_edges[0], e_edges[-1]),
            aspect="auto",
            cmap="viridis",
            vmin=-6.0,
            vmax=-2.0,
            interpolation="nearest",
        )
        if case.resonance_z is not None:
            ax.axvline(case.resonance_z, color="red", linestyle=":", linewidth=5.0)
        ax.axvline(8.0, color="limegreen", linewidth=5.0)

        x_m, b_m, _, _ = read_field_profile(case, case.snapshots[-1])
        ax2 = ax.twinx()
        ax2.plot(x_m, b_m, color="red", linewidth=4.0)
        ax2.set_ylabel("B0 [T]", color="red")
        ax2.tick_params(axis="y", colors="red", width=2.0, length=7)
        ax2.set_ylim(0.0, max(1.45, 1.05 * float(np.nanmax(b_m))))
        ax2.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

        hot_weight_fraction = (
            float(np.sum(phase.weights[hot_mask]) / np.sum(phase.weights))
            if np.sum(phase.weights) > 0.0
            else float("nan")
        )
        ax.set_xlim(-2.0, 8.0)
        ax.set_ylim(0.0, energy_max)
        ax.set_xlabel("z [m]")
        ax.set_ylabel("E_e [eV]")
        fig.suptitle(f"I: {current_a:.0f} [A] (PICOS++ nonrel)", y=0.985)
        ax.set_title(
            f"{case.label} | E_e >= {threshold_ev:.0f} eV | hot weight fraction {hot_weight_fraction:.2e}",
            pad=8,
        )
        ax.grid(True, color="0.7", alpha=0.3, linewidth=1.4)
        ax.tick_params(width=2.0, length=7)

        cax = inset_axes(
            ax,
            width="4.0%",
            height="76%",
            loc="upper right",
            bbox_to_anchor=(-0.13, -0.02, 1.0, 1.0),
            bbox_transform=ax.transAxes,
            borderpad=0,
        )
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label("log10(f_hot)", fontsize=22, labelpad=10)
        cbar.set_ticks([-6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2])
        cbar.ax.tick_params(width=2.0, length=7)

        fig.subplots_adjust(left=0.10, right=0.86, bottom=0.14, top=0.82)
        fig.savefig(out_path, dpi=220)
        plt.close(fig)


def plot_hot_filtered_fig17_maps(
    cases: Sequence[CaseInfo],
    phase_data: Dict[str, PhaseSpaceData],
    thresholds_ev: Sequence[float],
    energy_max: float,
    current_a: float,
    out_dir: Path,
) -> None:
    for case in cases:
        phase = phase_data[case.name]
        for threshold in thresholds_ev:
            hot_mask = phase.energy_ev >= threshold
            plot_single_fig17_hot_map(
                case,
                phase,
                hot_mask,
                energy_max,
                current_a,
                threshold,
                out_dir / f"fig17_hot_Egt{int(threshold)}_{case.name}.png",
            )
            trapped_mask = hot_mask & phase.trapped_magnetic
            plot_single_fig17_hot_map(
                case,
                phase,
                trapped_mask,
                energy_max,
                current_a,
                threshold,
                out_dir / f"fig17_magnetic_trapped_Egt{int(threshold)}_{case.name}.png",
            )


def plot_pitch_loss_cone(
    cases: Sequence[CaseInfo],
    phase_data: Dict[str, PhaseSpaceData],
    threshold_ev: float,
    out_dir: Path,
) -> None:
    fig, axes = plt.subplots(len(cases), 1, figsize=(12.0, 4.2 * len(cases)), sharex=True)
    if len(cases) == 1:
        axes = [axes]
    for ax, case in zip(axes, cases):
        phase = phase_data[case.name]
        hot = phase.energy_ev >= threshold_ev
        z_edges = np.linspace(-2.0, 8.0, 180)
        pitch_edges = np.linspace(0.0, 1.0, 120)
        hist, _, _ = np.histogram2d(
            phase.z[hot],
            phase.pitch_sin2[hot],
            bins=(z_edges, pitch_edges),
            weights=phase.weights[hot],
        )
        prob = hist / np.sum(hist) if np.sum(hist) > 0.0 else hist
        log_prob = np.full_like(prob, np.nan, dtype=np.float64)
        positive = prob > 0.0
        log_prob[positive] = np.log10(prob[positive])
        im = ax.imshow(
            np.ma.masked_invalid(log_prob.T),
            origin="lower",
            extent=(z_edges[0], z_edges[-1], pitch_edges[0], pitch_edges[-1]),
            aspect="auto",
            cmap="magma",
            vmin=-6.0,
            vmax=-2.0,
            interpolation="nearest",
        )
        x_m, b_m, _, _ = read_field_profile(case, case.snapshots[-1])
        b_left_max = np.maximum.accumulate(b_m)
        b_right_max = np.maximum.accumulate(b_m[::-1])[::-1]
        with np.errstate(divide="ignore", invalid="ignore"):
            loss = np.maximum(b_m / b_left_max, b_m / b_right_max)
        loss = np.clip(loss, 0.0, 1.0)
        ax.plot(x_m, loss, color="cyan", linewidth=2.5, label="two-sided magnetic loss cone")
        if case.resonance_z is not None:
            ax.axvline(case.resonance_z, color="white", linestyle="--", linewidth=1.5)
        ax.axvline(8.0, color="limegreen", linewidth=1.8)
        ax.set_ylabel("sin2(alpha)")
        ax.set_title(f"{case.label}: E_e >= {threshold_ev:.0f} eV")
        ax.legend(frameon=False, loc="lower right", fontsize=9)
        set_axis_style(ax)
    axes[-1].set_xlabel("z [m]")
    for ax in axes:
        ax.set_xlim(-2.0, 8.0)
        ax.set_ylim(0.0, 1.0)
    fig.subplots_adjust(left=0.08, right=0.84, bottom=0.08, top=0.93, hspace=0.30)
    cbar = fig.colorbar(im, ax=list(axes), pad=0.02, fraction=0.035)
    cbar.set_label("log10 weighted hot-population probability")
    fig.suptitle("Pitch angle versus magnetic loss cone")
    fig.savefig(out_dir / f"pitch_loss_cone_Egt{int(threshold_ev)}.png", dpi=220)
    plt.close(fig)


def plot_mirror_trapping_summary(
    cases: Sequence[CaseInfo],
    phase_data: Dict[str, PhaseSpaceData],
    threshold_ev: float,
    out_dir: Path,
) -> None:
    fig, axes = plt.subplots(len(cases), 2, figsize=(15.5, 4.0 * len(cases)), sharex=True)
    if len(cases) == 1:
        axes = np.asarray([axes])
    z_edges = np.linspace(-2.0, 8.0, 120)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    for row, case in enumerate(cases):
        phase = phase_data[case.name]
        hot = phase.energy_ev >= threshold_ev
        mag = hot & phase.trapped_magnetic
        eff = hot & phase.trapped_effective
        all_hot_density, _ = np.histogram(phase.z[hot], bins=z_edges, weights=phase.weights[hot])
        mag_density, _ = np.histogram(phase.z[mag], bins=z_edges, weights=phase.weights[mag])
        eff_density, _ = np.histogram(phase.z[eff], bins=z_edges, weights=phase.weights[eff])
        norm = np.nanmax(all_hot_density) if np.nanmax(all_hot_density) > 0.0 else 1.0
        ax = axes[row, 0]
        ax.plot(z_centers, all_hot_density / norm, color="0.25", linewidth=2.0, label="all hot")
        ax.plot(z_centers, mag_density / norm, color="#1f77b4", linewidth=2.0, label="magnetic trapped")
        ax.plot(z_centers, eff_density / norm, color="#d62728", linewidth=2.0, label="magnetic + phi trapped")
        ax.set_ylabel(case.label)
        ax.set_title(f"relative line density, E_e >= {threshold_ev:.0f} eV")
        ax.legend(frameon=False, fontsize=9)

        total_hot_by_bin, _ = np.histogram(phase.z[hot], bins=z_edges, weights=phase.weights[hot])
        mag_by_bin, _ = np.histogram(phase.z[mag], bins=z_edges, weights=phase.weights[mag])
        eff_by_bin, _ = np.histogram(phase.z[eff], bins=z_edges, weights=phase.weights[eff])
        with np.errstate(divide="ignore", invalid="ignore"):
            mag_frac = np.where(total_hot_by_bin > 0.0, mag_by_bin / total_hot_by_bin, np.nan)
            eff_frac = np.where(total_hot_by_bin > 0.0, eff_by_bin / total_hot_by_bin, np.nan)
        axes[row, 1].plot(z_centers, mag_frac, color="#1f77b4", linewidth=2.0, label="magnetic")
        axes[row, 1].plot(z_centers, eff_frac, color="#d62728", linewidth=2.0, label="magnetic + phi")
        axes[row, 1].set_ylim(-0.05, 1.05)
        axes[row, 1].set_title("hot-particle trapped fraction by z")
        axes[row, 1].legend(frameon=False, fontsize=9)

        for col in range(2):
            if case.resonance_z is not None:
                axes[row, col].axvline(case.resonance_z, color="0.45", linestyle="--", linewidth=1.3)
            axes[row, col].axvline(8.0, color="limegreen", linewidth=1.5)
            axes[row, col].set_xlim(-2.0, 8.0)
            set_axis_style(axes[row, col])
    axes[-1, 0].set_xlabel("z [m]")
    axes[-1, 1].set_xlabel("z [m]")
    fig.suptitle("Mirror-trapping diagnostic from final snapshot")
    fig.tight_layout()
    fig.savefig(out_dir / f"mirror_trapping_summary_Egt{int(threshold_ev)}.png", dpi=220)
    plt.close(fig)


def write_mirror_metrics(
    cases: Sequence[CaseInfo],
    phase_data: Dict[str, PhaseSpaceData],
    thresholds_ev: Sequence[float],
    out_dir: Path,
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for case in cases:
        phase = phase_data[case.name]
        total_weight = float(np.sum(phase.weights))
        for threshold in thresholds_ev:
            hot = phase.energy_ev >= threshold
            hot_weight = float(np.sum(phase.weights[hot]))
            if hot_weight > 0.0:
                mag_weight = float(np.sum(phase.weights[hot & phase.trapped_magnetic]))
                eff_weight = float(np.sum(phase.weights[hot & phase.trapped_effective]))
                pitch_mean = float(np.average(phase.pitch_sin2[hot], weights=phase.weights[hot]))
                loss_mean = float(np.average(phase.loss_cone_sin2[hot], weights=phase.weights[hot]))
                z_mean = float(np.average(phase.z[hot], weights=phase.weights[hot]))
            else:
                mag_weight = eff_weight = pitch_mean = loss_mean = z_mean = float("nan")
            rows.append(
                {
                    "case": case.name,
                    "label": case.label,
                    "threshold_eV": threshold,
                    "hot_particle_count": int(np.count_nonzero(hot)),
                    "hot_weight": hot_weight,
                    "hot_weight_fraction": hot_weight / total_weight if total_weight > 0.0 else float("nan"),
                    "magnetic_trapped_fraction_of_hot": mag_weight / hot_weight
                    if hot_weight > 0.0
                    else float("nan"),
                    "effective_trapped_fraction_of_hot": eff_weight / hot_weight
                    if hot_weight > 0.0
                    else float("nan"),
                    "weighted_mean_sin2_alpha_hot": pitch_mean,
                    "weighted_mean_loss_cone_sin2_hot": loss_mean,
                    "weighted_mean_z_hot_m": z_mean,
                }
            )
    write_rows(out_dir / "mirror_trapping_metrics.csv", rows)
    return rows


def plot_density_fields(
    cases: Sequence[CaseInfo],
    final_data: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]],
    out_dir: Path,
) -> None:
    fig, axes = plt.subplots(len(cases), 3, figsize=(15.5, 3.8 * len(cases)), sharex="col")
    if len(cases) == 1:
        axes = np.asarray([axes])
    z_edges = np.linspace(-2.0, 8.0, 220)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    for row, case in enumerate(cases):
        z, _, weights = final_data[case.name]
        hist, _ = np.histogram(z, bins=z_edges, weights=weights)
        dz = np.diff(z_edges)
        density = hist / dz
        density /= np.nanmax(density) if np.nanmax(density) > 0 else 1.0
        x_m, b_m, e_m, phi_m = read_field_profile(case, case.snapshots[-1])
        axes[row, 0].plot(z_centers, density, color=CASE_COLORS.get(case.name), linewidth=2.0)
        axes[row, 0].set_ylabel(case.label)
        axes[row, 0].set_title("relative electron line density")
        axes[row, 1].plot(x_m, e_m, color="0.2", linewidth=1.8)
        axes[row, 1].set_title("E_parallel")
        axes[row, 1].set_ylabel("E [V/m]")
        axes[row, 2].plot(x_m, phi_m, color="#9467bd", linewidth=1.8, label="phi")
        ax2 = axes[row, 2].twinx()
        ax2.plot(x_m, b_m, color="#d62728", linewidth=1.5, label="B")
        ax2.set_ylabel("B [T]", color="#d62728")
        ax2.tick_params(axis="y", colors="#d62728")
        axes[row, 2].set_title("phi and B")
        axes[row, 2].set_ylabel("phi [V]")
        for col in range(3):
            if case.resonance_z is not None:
                axes[row, col].axvline(case.resonance_z, color="0.5", linestyle="--", linewidth=1.2)
            axes[row, col].axvline(8.0, color="limegreen", linewidth=1.4)
            set_axis_style(axes[row, col])
    for col in range(3):
        axes[-1, col].set_xlabel("z [m]")
        axes[-1, col].set_xlim(-2.0, 8.0)
    fig.suptitle("Final particle density and field profiles")
    fig.tight_layout()
    fig.savefig(out_dir / "density_and_fields_final.png", dpi=220)
    plt.close(fig)


def plot_bfield_and_rf(cases: Sequence[CaseInfo], out_dir: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11.5, 8.0), sharex=True)
    for case in cases:
        x_m, b_m, _, _ = read_field_profile(case, case.snapshots[-1])
        axes[0].plot(x_m, b_m, label=case.label, linewidth=2.0, color=CASE_COLORS.get(case.name))
        if case.resonance_z is not None:
            axes[0].axvline(case.resonance_z, color=CASE_COLORS.get(case.name), linestyle="--", alpha=0.75)
        x1 = as_float(case.input_values, "RF_electron_x1")
        x2 = as_float(case.input_values, "RF_electron_x2")
        amp = as_float(case.input_values, "RF_electron_EfieldAmplitude") or 1.0
        if x1 is not None and x2 is not None:
            y = np.where((x_m >= x1) & (x_m <= x2), amp, 0.0)
            axes[1].plot(x_m, y, label=case.label, linewidth=2.0, color=CASE_COLORS.get(case.name))
    axes[0].set_ylabel("B [T]")
    axes[0].set_title("Magnetic field and estimated resonance markers")
    axes[1].set_xlabel("z [m]")
    axes[1].set_ylabel("RF E-field window [V/m]")
    axes[1].set_title("Configured electron RF window")
    for ax in axes:
        ax.axvline(8.0, color="limegreen", linewidth=1.4)
        ax.legend(frameon=False, fontsize=9)
        ax.set_xlim(-2.0, 8.0)
        set_axis_style(ax)
    fig.tight_layout()
    fig.savefig(out_dir / "bfield_and_rf_windows.png", dpi=220)
    plt.close(fig)


def read_profile_values(
    case: CaseInfo, input_key: str, fallback_suffix: Optional[str] = None
) -> Optional[np.ndarray]:
    if case.input_file is None:
        return None
    file_name = case.input_values.get(input_key)
    profile_path: Optional[Path] = None
    if file_name:
        profile_path = case.input_file.parent / file_name
    elif fallback_suffix:
        matches = sorted(case.input_file.parent.glob(f"*{fallback_suffix}"))
        profile_path = matches[0] if matches else None
    if profile_path is None:
        return None
    if not profile_path.exists():
        return None
    try:
        return np.asarray(np.loadtxt(profile_path), dtype=np.float64).reshape(-1)
    except Exception:
        return None


def plot_input_profiles(cases: Sequence[CaseInfo], out_dir: Path) -> None:
    if not cases:
        return
    case = cases[0]
    x_m, b_m, _, _ = read_field_profile(case, case.snapshots[0])
    profiles = [
        ("IC_ne_fileName", "normalized n_e", "_ne_norm.txt"),
        ("IC_Te_fileName", "normalized T_e", None),
        ("pairSource_fileName", "normalized pair source", None),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.4), sharex=True)
    axes = axes.reshape(-1)
    axes[0].plot(x_m, b_m, color="#d62728", linewidth=2.0)
    axes[0].set_title("B field used by PICOS++")
    axes[0].set_ylabel("B [T]")
    for ax_idx, (input_key, title, fallback_suffix) in enumerate(profiles, start=1):
        values = read_profile_values(case, input_key, fallback_suffix)
        if values is None:
            continue
        xx = x_m if len(values) == len(x_m) else np.linspace(-2.0, 8.0, len(values))
        axes[ax_idx].plot(xx, values, color="0.2", linewidth=2.0)
        axes[ax_idx].set_title(title)
        axes[ax_idx].set_ylabel("normalized value")

    for case_i in cases:
        for ax in axes:
            if case_i.resonance_z is not None:
                ax.axvline(
                    case_i.resonance_z,
                    color=CASE_COLORS.get(case_i.name),
                    linestyle="--",
                    linewidth=1.4,
                    alpha=0.65,
                )
            ax.axvline(8.0, color="limegreen", linewidth=1.2)
            ax.set_xlim(-2.0, 8.0)
            set_axis_style(ax)
    for ax in axes[-2:]:
        ax.set_xlabel("z [m]")
    fig.suptitle("Input profiles and resonance markers")
    fig.tight_layout()
    fig.savefig(out_dir / "input_profiles_and_resonances.png", dpi=220)
    plt.close(fig)


def write_case_inventory(cases: Sequence[CaseInfo], out_dir: Path) -> None:
    rows = []
    keys = [
        "SW_EfieldSolve",
        "SW_fieldSolveModel",
        "SW_Collisions",
        "SW_pairSource",
        "SW_RFheating",
        "SW_RFheatingElectrons",
        "SW_RFheatingIons",
        "RF_electron_Prf",
        "RF_electron_freq",
        "RF_electron_x1",
        "RF_electron_x2",
        "simulationTime",
        "outputCadence",
        "LX_min",
        "LX_max",
        "CV_ne",
        "CV_Te",
        "IC_ne",
        "IC_Te",
    ]
    for case in cases:
        row: Dict[str, object] = {
            "case": case.name,
            "label": case.label,
            "tag": case.tag,
            "hdf5_dir": str(case.hdf5_dir),
            "input_file": str(case.input_file) if case.input_file else "",
            "num_particle_files": len(case.particles),
            "num_snapshots": len(case.snapshots),
            "final_time_us": case.times_s[-1] * 1e6,
            "resonance_z_estimate_m": case.resonance_z if case.resonance_z is not None else "",
        }
        for key in keys:
            row[key] = case.input_values.get(key, "")
        rows.append(row)
    write_rows(out_dir / "case_inventory.csv", rows)


def write_readme(
    out_dir: Path,
    cases: Sequence[CaseInfo],
    final_rows: Sequence[Dict[str, object]],
) -> None:
    lines = [
        "# PICOS++ MPEX ECH Triplet Analysis",
        "",
        "This directory contains analysis products generated from the NERSC MPEX profile triplet archive.",
        "The EEDFs use `ions/spp_2` as the kinetic electron species and particle weights from `a_p`.",
        "",
        "## Cases",
    ]
    for case in cases:
        freq = case.rf_frequency_hz / 1e9 if case.rf_frequency_hz else float("nan")
        zres = case.resonance_z if case.resonance_z is not None else float("nan")
        lines.append(f"- {case.label}: f = {freq:.6g} GHz, resonance marker z = {zres:.6g} m")
    lines.extend(["", "## Final Snapshot Metrics"])
    for row in final_rows:
        lines.append(
            "- {label}: t = {time:.2f} us, mean = {mean:.2f} eV, p99 = {p99:.2f} eV, "
            "p99.9 = {p999:.2f} eV, max = {maxe:.2f} eV, weight fraction >300 eV = {f300:.3e}".format(
                label=row["label"],
                time=float(row["time_s"]) * 1e6,
                mean=float(row["mean_energy_eV"]),
                p99=float(row["p99_energy_eV"]),
                p999=float(row["p999_energy_eV"]),
                maxe=float(row["max_energy_eV"]),
                f300=float(row["frac_weight_E_gt_300eV"]),
            )
        )
    lines.extend(
        [
            "",
            "## Files",
            "- `eedf_final_global.png`: global final-time EEDF comparison plus initial EEDF.",
            "- `eedf_final_logbins_and_tail.png`: log-binned final EEDF and cumulative high-energy tail.",
            "- `eedf_time_evolution_*.png`: EEDF evolution for each case.",
            "- `eedf_final_regions.png`: final EEDF split by axial region.",
            "- `electron_energy_vs_z_final.png`: Fig. 17-style final energy-versus-z distributions with B-field overlay.",
            "- `fig17_style_*.png`: individual paper-style energy-versus-z maps with B-field overlay.",
            "- `fig17_hot_Egt*_*.png`: hot-electron-only Fig. 17-style maps.",
            "- `fig17_magnetic_trapped_Egt*_*.png`: Fig. 17-style maps filtered to magnetic-mirror-eligible electrons.",
            "- `pitch_loss_cone_Egt100.png`: hot-electron pitch-angle distribution versus magnetic loss cone.",
            "- `mirror_trapping_summary_Egt100.png`: hot-electron trapped fraction and line-density diagnostics.",
            "- `energy_metrics_vs_time.png`: weighted mean, p99, p99.9, and max energy versus time.",
            "- `tail_fraction_vs_time.png`: weighted high-energy tail fractions versus time.",
            "- `density_and_fields_final.png`: final relative line density, E_parallel, phi, and B profiles.",
            "- `bfield_and_rf_windows.png`: B-field and RF heating window setup.",
            "- `input_profiles_and_resonances.png`: normalized input profiles and resonance markers.",
            "- `time_metrics.csv`, `final_metrics.csv`, `case_inventory.csv`, `mirror_trapping_metrics.csv`: numeric summaries.",
        ]
    )
    (out_dir / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    root = args.root
    if root is None:
        if args.archive is None:
            raise SystemExit("Provide --root or --archive")
        root = extract_archive(args.archive)
    root = root.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    cases = discover_cases(root)
    write_case_inventory(cases, out_dir)

    print(f"Found {len(cases)} cases under {root}")
    for case in cases:
        print(
            f"  {case.name}: {len(case.particles)} particle files, "
            f"{len(case.snapshots)} snapshots, final t={case.times_s[-1] * 1e6:.3f} us"
        )

    time_rows: List[Dict[str, object]] = []
    final_data: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    final_rows: List[Dict[str, object]] = []

    for case in cases:
        for snap_index, snap in enumerate(case.snapshots):
            z, energy, weights = read_particles(case, snap, args.electron_species)
            row = summarize_snapshot(case, snap_index, z, energy, weights)
            time_rows.append(row)
            if snap_index == len(case.snapshots) - 1:
                final_data[case.name] = (z, energy, weights)
                final_rows.append(row)
        print(f"  analyzed {case.name}")

    write_rows(out_dir / "time_metrics.csv", time_rows)
    write_rows(out_dir / "final_metrics.csv", final_rows)

    if args.energy_max is None:
        global_max = max(float(np.max(final_data[case.name][1])) for case in cases)
        energy_max = max(global_max * 1.05, 1000.0)
        energy_max = math.ceil(energy_max / 100.0) * 100.0
    else:
        energy_max = args.energy_max
    bins = np.linspace(0.0, energy_max, args.energy_bins + 1)
    hot_thresholds = [float(item) for item in args.hot_thresholds.split(",") if item.strip()]
    mirror_thresholds = sorted(set([0.0, 100.0, 300.0, 500.0, 1000.0] + hot_thresholds))

    phase_data = {
        case.name: read_phase_space(case, case.snapshots[-1], args.electron_species) for case in cases
    }
    mirror_rows = write_mirror_metrics(cases, phase_data, mirror_thresholds, out_dir)

    plot_bfield_and_rf(cases, out_dir)
    plot_input_profiles(cases, out_dir)
    plot_final_eedf(cases, final_data, bins, args.electron_species, out_dir)
    plot_final_eedf_logbins(cases, final_data, energy_max, args.electron_species, out_dir)
    plot_time_evolution_eedf(cases, args.electron_species, bins, out_dir)
    plot_region_eedf(cases, final_data, bins, out_dir)
    plot_energy_metrics(time_rows, out_dir)
    plot_tail_fraction(time_rows, out_dir)
    plot_energy_vs_z(cases, final_data, energy_max, args.z_bins, out_dir)
    plot_fig17_style_maps(cases, final_data, args.fig17_energy_max, args.figure_current_a, out_dir)
    plot_hot_filtered_fig17_maps(
        cases, phase_data, hot_thresholds, args.fig17_energy_max, args.figure_current_a, out_dir
    )
    if hot_thresholds:
        diagnostic_threshold = hot_thresholds[0]
        plot_pitch_loss_cone(cases, phase_data, diagnostic_threshold, out_dir)
        plot_mirror_trapping_summary(cases, phase_data, diagnostic_threshold, out_dir)
    plot_density_fields(cases, final_data, out_dir)
    write_readme(out_dir, cases, final_rows)

    print(f"Wrote analysis outputs to {out_dir}")


if __name__ == "__main__":
    main()
