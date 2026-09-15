#!/usr/bin/env python3
"""Analyze a two-stage MPEX steady-state then ECH PICOS++ restart run.

The NERSC restart bundles used for the MPEX ECH tests write HDF5 directly under
``picosFILES/outputFiles/HDF5``.  Some older PICOS analysis scripts assume a
named run directory below ``outputFiles``; this script keeps the two stages
separate, checks restart timing, and produces profile/EEDF diagnostics.
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


E_CHARGE = 1.602176634e-19
M_E = 9.1093837015e-31

SPECIES_LABEL = {
    "spp_1": "D+",
    "spp_2": "e-",
}


@dataclass
class PicosStage:
    key: str
    label: str
    hdf5_dir: Path
    input_dir: Path | None
    log_files: list[Path]
    params: dict[str, str]
    z: np.ndarray
    times: np.ndarray
    fields: dict[str, np.ndarray]
    moments: dict[str, dict[str, np.ndarray]]
    species: list[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path, help="Extracted root or .tar.gz archive.")
    parser.add_argument("--extract-dir", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=Path("validation/mpex_ech_restart_analysis"))
    parser.add_argument("--electron-species", default="spp_2")
    parser.add_argument("--ion-species", default="spp_1")
    parser.add_argument("--energy-max", type=float, default=3500.0)
    parser.add_argument("--smooth-window", type=int, default=3)
    return parser.parse_args()


def maybe_extract(path: Path, extract_dir: Path | None) -> Path:
    path = path.expanduser()
    if path.is_file() and path.suffixes[-2:] == [".tar", ".gz"]:
        if extract_dir is None:
            extract_dir = Path("/private/tmp") / f"{path.stem.replace('.tar', '')}_extract"
        extract_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(path, "r:gz") as tf:
            tf.extractall(extract_dir)
        return extract_dir
    return path


def numeric_keys(handle: h5py.File) -> list[str]:
    return sorted([key for key in handle.keys() if key.isdigit()], key=lambda item: int(item))


def read_1d(handle: h5py.File, path: str) -> np.ndarray | None:
    if path not in handle:
        return None
    return np.asarray(handle[path], dtype=np.float64).reshape(-1)


def read_scalar(handle: h5py.File, path: str, default: float = math.nan) -> float:
    values = read_1d(handle, path)
    if values is None or values.size == 0:
        return default
    return float(values[0])


def parse_key_values(path: Path | None) -> dict[str, str]:
    if path is None or not path.exists():
        return {}
    values: dict[str, str] = {}
    for raw in path.read_text(errors="ignore").splitlines():
        line = raw.split("//", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 2:
            values[parts[0]] = parts[1]
    return values


def find_input_dir(hdf5_dir: Path) -> Path | None:
    for parent in hdf5_dir.parents:
        candidate = parent / "picosFILES" / "inputFiles"
        if candidate.is_dir():
            return candidate
        if parent.name == "picosFILES":
            candidate = parent / "inputFiles"
            if candidate.is_dir():
                return candidate
    return None


def find_log_files(hdf5_dir: Path) -> list[Path]:
    for parent in hdf5_dir.parents:
        candidate = parent / "logs"
        if candidate.is_dir():
            return sorted(candidate.glob("*.log"))
    return []


def stage_key_from_path(hdf5_dir: Path) -> tuple[str, str]:
    text = str(hdf5_dir)
    if "steady_interactive" in text or "steady_2ms" in text:
        return "steady", "RF off steady stage"
    if "well_min" in text or "ech" in text.lower():
        return "ech", "well-minimum ECH restart"
    return hdf5_dir.parent.name, hdf5_dir.parent.name


def discover_hdf5_dirs(root: Path) -> dict[str, Path]:
    dirs = sorted({candidate.parent.resolve() for candidate in root.rglob("HDF5/main.h5")})
    if not dirs and (root / "main.h5").is_file():
        dirs = [root.resolve()]
    if not dirs:
        raise FileNotFoundError(f"No HDF5/main.h5 found under {root}")
    out: dict[str, Path] = {}
    for hdf5_dir in dirs:
        key, _ = stage_key_from_path(hdf5_dir)
        out[key] = hdf5_dir
    return out


def load_stage(hdf5_dir: Path) -> PicosStage:
    key, label = stage_key_from_path(hdf5_dir)
    input_dir = find_input_dir(hdf5_dir)
    input_file = input_dir / "input_file.input" if input_dir else None
    params = parse_key_values(input_file)

    with h5py.File(hdf5_dir / "main.h5", "r") as handle:
        z = read_1d(handle, "geometry/x_m")
        if z is None:
            z = read_1d(handle, "geometry/xAxis")
        if z is None:
            raise KeyError(f"{hdf5_dir / 'main.h5'} has no mesh coordinate")

    fields: dict[str, np.ndarray] = {}
    f0 = hdf5_dir / "FIELDS_FILE_0.h5"
    with h5py.File(f0, "r") as handle:
        steps = numeric_keys(handle)
        times = np.asarray([read_scalar(handle, f"{step}/time") for step in steps], dtype=np.float64)
        for name in ("BX_m", "dBX_m", "ddBX_m", "EX_m", "Phi_m"):
            cols: list[np.ndarray] = []
            for step in steps:
                values = read_1d(handle, f"{step}/fields/{name}/x")
                if values is None:
                    cols = []
                    break
                cols.append(values)
            if cols:
                fields[name] = np.column_stack(cols)

    moments: dict[str, dict[str, np.ndarray]] = {}
    p0 = hdf5_dir / "PARTICLES_FILE_0.h5"
    with h5py.File(p0, "r") as handle:
        psteps = numeric_keys(handle)
        if psteps != steps:
            steps = psteps
            times = np.asarray([read_scalar(handle, f"{step}/time") for step in steps], dtype=np.float64)
        species = sorted(handle[f"{steps[0]}/ions"].keys(), key=lambda item: int(item.split("_")[-1]))
        for spp in species:
            moments[spp] = {}
            for name in ("n_m", "Tpar_m", "Tper_m"):
                cols = []
                for step in steps:
                    values = read_1d(handle, f"{step}/ions/{spp}/{name}")
                    if values is None:
                        cols = []
                        break
                    cols.append(values)
                if cols:
                    moments[spp][name] = np.column_stack(cols)
            cols = []
            for step in steps:
                values = read_1d(handle, f"{step}/ions/{spp}/u_m/x")
                if values is None:
                    cols = []
                    break
                cols.append(values)
            if cols:
                moments[spp]["u_m"] = np.column_stack(cols)

    n = min(
        [z.size]
        + [arr.shape[0] for group in moments.values() for arr in group.values()]
        + [arr.shape[0] for arr in fields.values()]
    )
    return PicosStage(
        key=key,
        label=label,
        hdf5_dir=hdf5_dir,
        input_dir=input_dir,
        log_files=find_log_files(hdf5_dir),
        params=params,
        z=z[:n],
        times=times,
        fields={k: v[:n, :] for k, v in fields.items()},
        moments={spp: {k: v[:n, :] for k, v in group.items()} for spp, group in moments.items()},
        species=species,
    )


def moving_average(values: np.ndarray, window: int, axis: int = 0) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if window <= 1:
        return arr
    kernel = np.ones(window, dtype=np.float64) / float(window)
    return np.apply_along_axis(lambda line: np.convolve(line, kernel, mode="same"), axis, arr)


def mask_by_density(values: np.ndarray, density: np.ndarray, floor_fraction: float = 1.0e-3) -> np.ndarray:
    out = np.asarray(values, dtype=np.float64).copy()
    den = np.asarray(density, dtype=np.float64)
    peak = float(np.nanmax(den)) if np.any(np.isfinite(den)) else 0.0
    if peak > 0.0:
        out[den < peak * floor_fraction] = np.nan
    return out


def finite_window_mean(values: np.ndarray, indices: Iterable[int]) -> np.ndarray:
    idx = list(indices)
    if not idx:
        idx = [values.shape[1] - 1]
    subset = np.asarray(values[:, idx], dtype=np.float64)
    return np.nanmean(subset, axis=1)


def density_weighted_mean(z: np.ndarray, values: np.ndarray, density: np.ndarray) -> float:
    mask = np.isfinite(z) & np.isfinite(values) & np.isfinite(density) & (density > 0.0)
    if np.count_nonzero(mask) < 2:
        return math.nan
    numerator = np.trapz(values[mask] * density[mask], z[mask])
    denominator = np.trapz(density[mask], z[mask])
    return float(numerator / max(denominator, 1.0e-300))


def line_integral(z: np.ndarray, values: np.ndarray) -> float:
    mask = np.isfinite(z) & np.isfinite(values)
    if np.count_nonzero(mask) < 2:
        return math.nan
    return float(np.trapz(values[mask], z[mask]))


def rel_l2(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if not np.any(mask):
        return math.nan
    return float(np.linalg.norm(a[mask] - b[mask]) / max(np.linalg.norm(b[mask]), 1.0e-300))


def load_profile(input_dir: Path | None, file_name: str | None) -> np.ndarray | None:
    if input_dir is None or not file_name or file_name == "none":
        return None
    path = input_dir / file_name
    if not path.exists():
        return None
    return np.asarray(np.loadtxt(path), dtype=np.float64).reshape(-1)


def profile_axis(params: dict[str, str], n: int) -> np.ndarray:
    xmin = float(params.get("LX_min", "0.0"))
    xmax = float(params.get("LX_max", str(n - 1)))
    return np.linspace(xmin, xmax, n)


def interp_input(stage: PicosStage, key: str) -> np.ndarray | None:
    profile = load_profile(stage.input_dir, stage.params.get(key))
    if profile is None:
        return None
    xp = profile_axis(stage.params, profile.size)
    return np.interp(stage.z, xp, profile)


def summary_rows(stage: PicosStage, final_window: int = 5) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    final_start = max(0, stage.times.size - final_window)
    previous_start = max(0, final_start - final_window)
    final_idx = range(final_start, stage.times.size)
    previous_idx = range(previous_start, final_start)
    for spp in stage.species:
        label = SPECIES_LABEL.get(spp, spp)
        for var in ("n_m", "Tpar_m", "Tper_m", "u_m"):
            values = stage.moments.get(spp, {}).get(var)
            if values is None:
                continue
            initial = values[:, 0]
            final = finite_window_mean(values, final_idx)
            previous = finite_window_mean(values, previous_idx) if final_start > previous_start else initial
            peak_idx = int(np.nanargmax(final)) if np.any(np.isfinite(final)) else 0
            rows.append(
                {
                    "stage": stage.key,
                    "species": label,
                    "variable": var,
                    "snapshots": stage.times.size,
                    "t_first_ms": stage.times[0] * 1.0e3,
                    "t_last_ms": stage.times[-1] * 1.0e3,
                    "line_initial": line_integral(stage.z, initial),
                    "line_final_window": line_integral(stage.z, final),
                    "final_over_initial": line_integral(stage.z, final) / max(line_integral(stage.z, initial), 1.0e-300),
                    "late_window_l2_change": rel_l2(final, previous),
                    "last_step_l2_change": rel_l2(values[:, -1], values[:, -2]) if values.shape[1] > 1 else math.nan,
                    "final_peak": float(np.nanmax(final)) if np.any(np.isfinite(final)) else math.nan,
                    "final_peak_z_m": float(stage.z[peak_idx]) if stage.z.size else math.nan,
                }
            )
    return rows


def active_region_rows(stage: PicosStage) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    samples = [("start", 0), ("final", stage.times.size - 1)]
    for spp in stage.species:
        label = SPECIES_LABEL.get(spp, spp)
        group = stage.moments.get(spp, {})
        density = group.get("n_m")
        if density is None:
            continue
        for name, idx in samples:
            n = density[:, idx]
            active = n > max(float(np.nanmax(n)), 0.0) * 1.0e-3
            row: dict[str, object] = {
                "stage": stage.key,
                "sample": name,
                "species": label,
                "time_ms": stage.times[idx] * 1.0e3,
                "line_density_m-2": line_integral(stage.z, n),
                "peak_density_m-3": float(np.nanmax(n)) if np.any(np.isfinite(n)) else math.nan,
                "active_z_min_m": float(np.nanmin(stage.z[active])) if np.any(active) else math.nan,
                "active_z_max_m": float(np.nanmax(stage.z[active])) if np.any(active) else math.nan,
            }
            for var in ("Tpar_m", "Tper_m", "u_m"):
                values = group.get(var)
                if values is not None:
                    row[f"density_weighted_{var}"] = density_weighted_mean(stage.z, values[:, idx], n)
            rows.append(row)
    return rows


def read_particle_phase(stage: PicosStage, species: str, step: str) -> dict[str, np.ndarray]:
    z_chunks: list[np.ndarray] = []
    energy_chunks: list[np.ndarray] = []
    weight_chunks: list[np.ndarray] = []
    vpar_chunks: list[np.ndarray] = []
    vperp_chunks: list[np.ndarray] = []
    b_chunks: list[np.ndarray] = []
    base = f"{step}/ions/{species}"
    for path in sorted(stage.hdf5_dir.glob("PARTICLES_FILE_*.h5"), key=lambda p: int(p.stem.rsplit("_", 1)[1])):
        with h5py.File(path, "r") as handle:
            if base not in handle:
                continue
            group = handle[base]
            x = np.asarray(group["X_p"], dtype=np.float64).reshape(-1)
            weight = np.asarray(group["a_p"], dtype=np.float64).reshape(-1)
            velocity = np.asarray(group["V_p"], dtype=np.float64)
            if velocity.shape[0] != 2:
                velocity = velocity.T
            vpar = velocity[0]
            vperp = velocity[1]
            mass = M_E if species == "spp_2" else 3.321076e-27
            energy = 0.5 * mass * (vpar * vpar + vperp * vperp) / E_CHARGE
            b = np.asarray(group.get("BX_p", np.full_like(x, np.nan)), dtype=np.float64).reshape(-1)
            valid = np.isfinite(x) & np.isfinite(energy) & np.isfinite(weight) & (weight > 0.0)
            z_chunks.append(x[valid])
            energy_chunks.append(energy[valid])
            weight_chunks.append(weight[valid])
            vpar_chunks.append(vpar[valid])
            vperp_chunks.append(vperp[valid])
            b_chunks.append(b[valid])
    if not z_chunks:
        raise KeyError(f"No particle state found for {species} at snapshot {step}")
    return {
        "z": np.concatenate(z_chunks),
        "energy_eV": np.concatenate(energy_chunks),
        "weight": np.concatenate(weight_chunks),
        "vpar": np.concatenate(vpar_chunks),
        "vperp": np.concatenate(vperp_chunks),
        "B": np.concatenate(b_chunks),
    }


def weighted_percentile(values: np.ndarray, weights: np.ndarray, percentile: float) -> float:
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(mask):
        return math.nan
    x = values[mask]
    w = weights[mask]
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    cdf = np.cumsum(w)
    return float(np.interp(percentile / 100.0 * cdf[-1], cdf, x))


def particle_summary(stage: PicosStage, species: str, step: str, name: str) -> dict[str, object]:
    data = read_particle_phase(stage, species, step)
    energy = data["energy_eV"]
    weight = data["weight"]
    z = data["z"]
    total_w = float(np.sum(weight))
    return {
        "stage": stage.key,
        "snapshot_name": name,
        "snapshot": step,
        "time_ms": float(stage.times[int(step)] * 1.0e3),
        "macro_particles": int(energy.size),
        "sum_weight": total_w,
        "mean_energy_eV": float(np.average(energy, weights=weight)) if total_w > 0 else math.nan,
        "p50_energy_eV": weighted_percentile(energy, weight, 50.0),
        "p90_energy_eV": weighted_percentile(energy, weight, 90.0),
        "p99_energy_eV": weighted_percentile(energy, weight, 99.0),
        "p999_energy_eV": weighted_percentile(energy, weight, 99.9),
        "frac_weight_E_gt_100eV": float(np.sum(weight[energy >= 100.0]) / total_w) if total_w > 0 else math.nan,
        "frac_weight_E_gt_300eV": float(np.sum(weight[energy >= 300.0]) / total_w) if total_w > 0 else math.nan,
        "frac_weight_E_gt_1000eV": float(np.sum(weight[energy >= 1000.0]) / total_w) if total_w > 0 else math.nan,
        "weighted_z_mean_m": float(np.average(z, weights=weight)) if total_w > 0 else math.nan,
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_profiles(stages: dict[str, PicosStage], out_dir: Path, smooth_window: int) -> None:
    steady = stages.get("steady")
    ech = stages.get("ech")
    if steady is None or ech is None:
        return
    z = steady.z
    fig, axes = plt.subplots(3, 2, figsize=(16.5, 13.0), sharex=True)
    axes = axes.reshape(3, 2)

    for col, spp in enumerate(("spp_1", "spp_2")):
        label = SPECIES_LABEL.get(spp, spp)
        ax = axes[0, col]
        for stage, color, ls, name in (
            (steady, "0.4", "-", "steady start"),
            (steady, "k", "--", "steady final"),
            (ech, "#1f77b4", "-", "ECH start"),
            (ech, "#d62728", "-", "ECH final"),
        ):
            values = stage.moments[spp]["n_m"][:, 0 if "start" in name else -1]
            ax.plot(stage.z, moving_average(values, smooth_window), color=color, ls=ls, lw=2.0, label=name)
        source = interp_input(ech, "pairSource_fileName")
        if source is not None and col == 0:
            ymax = ax.get_ylim()[1]
            if np.nanmax(source) > 0 and ymax > 0:
                ax.plot(ech.z, source / np.nanmax(source) * ymax * 0.9, color="#2ca02c", lw=1.8, ls=":", label="source shape")
        ax.set_title(f"{label} density")
        ax.set_ylabel("n [m$^{-3}$]")
        ax.grid(True, alpha=0.25)
        ax.legend(frameon=False, fontsize=9)

        ax = axes[1, col]
        for stage, color, ls, name in (
            (steady, "k", "--", "steady final"),
            (ech, "#1f77b4", "-", "ECH start"),
            (ech, "#d62728", "-", "ECH final"),
        ):
            idx = 0 if name.endswith("start") else -1
            for var, alpha in (("Tpar_m", 1.0), ("Tper_m", 0.65)):
                values = stage.moments[spp][var][:, idx]
                values = mask_by_density(values, stage.moments[spp]["n_m"][:, idx])
                ax.plot(
                    stage.z,
                    moving_average(values, smooth_window),
                    color=color,
                    ls=ls if var == "Tpar_m" else ":",
                    lw=2.0,
                    alpha=alpha,
                    label=f"{name} {'Tpar' if var == 'Tpar_m' else 'Tper'}",
                )
        ax.set_title(f"{label} temperature")
        ax.set_ylabel("T [eV]")
        ax.grid(True, alpha=0.25)
        ax.legend(frameon=False, fontsize=8, ncol=2)

        ax = axes[2, col]
        for stage, color, ls, name in (
            (steady, "k", "--", "steady final"),
            (ech, "#1f77b4", "-", "ECH start"),
            (ech, "#d62728", "-", "ECH final"),
        ):
            values = stage.moments[spp]["u_m"][:, 0 if name.endswith("start") else -1]
            ax.plot(stage.z, moving_average(values, smooth_window), color=color, ls=ls, lw=2.0, label=name)
        ax.set_title(f"{label} parallel flow")
        ax.set_ylabel("u$_\\parallel$ [m/s]")
        ax.set_xlabel("z [m]")
        ax.grid(True, alpha=0.25)
        ax.legend(frameon=False, fontsize=9)

    for ax in axes.reshape(-1):
        ax.set_xlim(-2.0, 8.0)
    fig.suptitle("MPEX PICOS++ steady restart and 100 us ECH profile comparison")
    fig.tight_layout()
    fig.savefig(out_dir / "mpex_ech0914_profile_comparison.png", dpi=220)
    plt.close(fig)


def plot_density_zoom(stages: dict[str, PicosStage], out_dir: Path, smooth_window: int) -> None:
    steady = stages.get("steady")
    ech = stages.get("ech")
    if steady is None or ech is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(15.2, 5.3), sharex=True)
    for ax, spp in zip(axes, ("spp_1", "spp_2")):
        label = SPECIES_LABEL.get(spp, spp)
        curves = [
            (steady.z, steady.moments[spp]["n_m"][:, -1], "k", "--", "steady final / restart"),
            (ech.z, ech.moments[spp]["n_m"][:, 0], "#1f77b4", "-", "ECH start"),
            (ech.z, ech.moments[spp]["n_m"][:, -1], "#d62728", "-", "ECH final"),
        ]
        for z, values, color, ls, name in curves:
            ax.plot(z, moving_average(values, smooth_window), color=color, ls=ls, lw=2.3, label=name)
        ax.set_title(f"{label} density zoom, excluding t=0 IC")
        ax.set_ylabel("n [m$^{-3}$]")
        ax.set_xlabel("z [m]")
        ax.set_xlim(-2.0, 8.0)
        ax.grid(True, alpha=0.25)
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_dir / "mpex_ech0914_density_zoom_no_initial.png", dpi=220)
    plt.close(fig)


def plot_time_maps(stage: PicosStage, out_dir: Path, smooth_window: int, skip_initial: bool) -> None:
    panels: list[tuple[str, np.ndarray, str, str]] = []
    for spp in ("spp_1", "spp_2"):
        if spp not in stage.moments:
            continue
        label = SPECIES_LABEL.get(spp, spp)
        density = stage.moments[spp]["n_m"]
        panels.append((f"{label} density", density, "viridis", "n [m$^{-3}$]"))
        temp = stage.moments[spp]["Tpar_m"].copy()
        for jj in range(temp.shape[1]):
            temp[:, jj] = mask_by_density(temp[:, jj], density[:, jj])
        panels.append((f"{label} Tpar", temp, "magma", "Tpar [eV]"))
    if not panels:
        return
    fig, axes = plt.subplots(2, 2, figsize=(15.8, 10.0), sharex=True, squeeze=False)
    tidx = np.arange(stage.times.size)
    suffix = "no_initial" if skip_initial and stage.times.size > 1 else "with_initial"
    if skip_initial and stage.times.size > 1:
        tidx = tidx[1:]
    t_ms = stage.times[tidx] * 1.0e3
    for ax, (title, values, cmap, label) in zip(axes.reshape(-1), panels):
        smoothed = moving_average(values[:, tidx], smooth_window)
        finite = smoothed[np.isfinite(smoothed)]
        vmax = float(np.nanpercentile(finite, 99.5)) if finite.size else 1.0
        if not np.isfinite(vmax) or vmax <= 0:
            vmax = 1.0
        im = ax.pcolormesh(stage.z, t_ms, smoothed.T, shading="auto", cmap=cmap, vmin=0.0, vmax=vmax)
        ax.set_title(f"{stage.label}: {title}")
        ax.set_ylabel("time [ms]")
        ax.grid(True, alpha=0.18)
        cbar = fig.colorbar(im, ax=ax, pad=0.01)
        cbar.set_label(label)
    for ax in axes[-1, :]:
        ax.set_xlabel("z [m]")
    for ax in axes.reshape(-1):
        ax.set_xlim(-2.0, 8.0)
    fig.tight_layout()
    fig.savefig(out_dir / f"mpex_ech0914_{stage.key}_time_maps_{suffix}.png", dpi=220)
    plt.close(fig)


def plot_fields(stages: dict[str, PicosStage], out_dir: Path, smooth_window: int) -> None:
    ech = stages.get("ech")
    steady = stages.get("steady")
    if ech is None:
        return
    fig, axes = plt.subplots(2, 1, figsize=(13.5, 9.0), sharex=True)
    if "BX_m" in ech.fields:
        axes[0].plot(ech.z, moving_average(ech.fields["BX_m"][:, -1], smooth_window), "k", lw=2.6, label="B")
        try:
            x1 = float(ech.params.get("RF_electron_x1", "nan"))
            axes[0].axvline(x1 + 0.02, color="red", lw=2.2, ls=":", label="ECH resonance")
        except ValueError:
            pass
        axes[0].set_ylabel("B [T]")
        axes[0].grid(True, alpha=0.25)
        axes[0].legend(frameon=False)
    for stage, color, name in ((steady, "0.45", "steady final"), (ech, "#d62728", "ECH final")):
        if stage is None or "EX_m" not in stage.fields:
            continue
        axes[1].plot(stage.z, moving_average(stage.fields["EX_m"][:, -1], smooth_window), color=color, lw=2.2, label=name)
    axes[1].set_ylabel("E$_\\parallel$ [V/m]")
    axes[1].set_xlabel("z [m]")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(frameon=False)
    axes[1].set_xlim(-2.0, 8.0)
    fig.suptitle("MPEX PICOS++ magnetic and electrostatic fields")
    fig.tight_layout()
    fig.savefig(out_dir / "mpex_ech0914_fields.png", dpi=220)
    plt.close(fig)


def eedf_histogram(energy: np.ndarray, weights: np.ndarray, bins: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    hist, edges = np.histogram(energy, bins=bins, weights=weights)
    widths = np.diff(edges)
    total = np.sum(hist)
    pdf = hist / max(total, 1.0e-300) / widths
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, pdf


def plot_eedf(stages: dict[str, PicosStage], out_dir: Path, species: str, energy_max: float) -> None:
    steady = stages.get("steady")
    ech = stages.get("ech")
    if steady is None or ech is None:
        return
    samples = [
        ("steady start", steady, "0", "0.45"),
        ("steady final / ECH start", steady, str(steady.times.size - 1), "k"),
        ("ECH final", ech, str(ech.times.size - 1), "#d62728"),
    ]
    bins = np.linspace(0.0, energy_max, 240)
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.4), sharex=True)
    particle_rows: list[dict[str, object]] = []
    for label, stage, step, color in samples:
        data = read_particle_phase(stage, species, step)
        centers, pdf = eedf_histogram(data["energy_eV"], data["weight"], bins)
        axes[0].plot(centers, pdf, lw=2.3, color=color, label=label)
        axes[1].semilogy(centers, np.maximum(pdf, 1.0e-30), lw=2.3, color=color, label=label)
        particle_rows.append(particle_summary(stage, species, step, label))
    for ax in axes:
        ax.set_xlabel("electron kinetic energy [eV]")
        ax.grid(True, alpha=0.25)
        ax.legend(frameon=False)
    axes[0].set_ylabel("weighted PDF [1/eV]")
    axes[0].set_title("linear EEDF")
    axes[1].set_ylabel("weighted PDF [1/eV]")
    axes[1].set_title("semilog EEDF")
    fig.tight_layout()
    fig.savefig(out_dir / "mpex_ech0914_eedf.png", dpi=220)
    plt.close(fig)
    write_csv(out_dir / "mpex_ech0914_particle_energy_summary.csv", particle_rows)


def plot_energy_z(stage: PicosStage, out_dir: Path, species: str, energy_max: float) -> None:
    step = str(stage.times.size - 1)
    data = read_particle_phase(stage, species, step)
    z_edges = np.linspace(-2.0, 8.0, 220)
    e_edges = np.linspace(0.0, energy_max, 260)
    hist, _, _ = np.histogram2d(data["z"], data["energy_eV"], bins=(z_edges, e_edges), weights=data["weight"])
    density = hist.T
    positive = density[density > 0]
    log_density = np.full_like(density, np.nan, dtype=np.float64)
    if positive.size:
        log_density[density > 0] = np.log10(density[density > 0] / np.nanmax(positive))

    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    e_centers = 0.5 * (e_edges[:-1] + e_edges[1:])
    fig, axes = plt.subplots(2, 1, figsize=(14.5, 10.0), sharex=True)
    im0 = axes[0].pcolormesh(z_centers, e_centers, np.nan_to_num(density, nan=0.0), shading="auto", cmap="viridis")
    axes[0].set_title(f"{stage.label}: weighted electron energy map, linear")
    axes[0].set_ylabel("E$_e$ [eV]")
    cbar0 = fig.colorbar(im0, ax=axes[0], pad=0.01)
    cbar0.set_label("weighted counts")

    im1 = axes[1].pcolormesh(z_centers, e_centers, log_density, shading="auto", cmap="viridis", vmin=-6.0, vmax=0.0)
    axes[1].set_title("Fig. 17-style log$_{10}$ normalized f(E,z)")
    axes[1].set_ylabel("E$_e$ [eV]")
    axes[1].set_xlabel("z [m]")
    cbar1 = fig.colorbar(im1, ax=axes[1], pad=0.01)
    cbar1.set_label("log$_{10}$(f/f$_{max}$)")

    if "BX_m" in stage.fields:
        b = moving_average(stage.fields["BX_m"][:, -1], 3)
        b_scale = energy_max / max(np.nanmax(b), 1.0e-300)
        for ax in axes:
            ax.plot(stage.z, b * b_scale, color="red", lw=2.4, label="B scaled")
            try:
                ax.axvline(float(stage.params.get("RF_electron_x1", "nan")) + 0.02, color="red", ls=":", lw=2.0)
            except ValueError:
                pass
            ax.grid(True, alpha=0.2)
            ax.legend(frameon=False)
            ax.set_xlim(-2.0, 8.0)
            ax.set_ylim(0.0, energy_max)
    fig.tight_layout()
    fig.savefig(out_dir / f"mpex_ech0914_{stage.key}_energy_vs_z.png", dpi=220)
    plt.close(fig)


def write_summary_markdown(
    out_dir: Path,
    stages: dict[str, PicosStage],
    moment_rows: list[dict[str, object]],
    active_rows: list[dict[str, object]],
    particle_rows: list[dict[str, object]],
) -> None:
    steady = stages.get("steady")
    ech = stages.get("ech")
    lines: list[str] = ["# MPEX ECH Restart Analysis", ""]
    if steady:
        lines.append(
            f"- Steady snapshots: {steady.times.size}, saved time range {steady.times[0] * 1e3:.6f} to {steady.times[-1] * 1e3:.6f} ms."
        )
    if ech:
        lines.append(
            f"- ECH snapshots: {ech.times.size}, saved time range {ech.times[0] * 1e3:.6f} to {ech.times[-1] * 1e3:.6f} ms."
        )
        lines.append(
            f"- ECH flags: SW_RFheating={ech.params.get('SW_RFheating')}, SW_RFheatingElectrons={ech.params.get('SW_RFheatingElectrons')}, SW_RFheatingIons={ech.params.get('SW_RFheatingIons')}, restart_enabled={ech.params.get('restart_enabled')}."
        )
        lines.append(
            f"- Pair source: rate={ech.params.get('pairSource_rate')} m^-3 s^-1, mean z={ech.params.get('pairSource_mean_x')} m, sigma={ech.params.get('pairSource_sigma_x')} m."
        )
    if steady and ech:
        lines.append(
            f"- Restart continuity: ECH first saved time minus steady last saved time = {(ech.times[0] - steady.times[-1]) * 1e6:.3f} us."
        )
    lines.append("")
    lines.append("## Density Moment Metrics")
    lines.append("")
    lines.append("| stage | species | final/initial line density | late-window L2 | last-step L2 | final peak z [m] | final peak density [m^-3] |")
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for row in moment_rows:
        if row["variable"] != "n_m" or row["species"] not in ("D+", "e-"):
            continue
        lines.append(
            "| {stage} | {species} | {ratio:.4g} | {late:.4g} | {last:.4g} | {z:.4g} | {peak:.4e} |".format(
                stage=row["stage"],
                species=row["species"],
                ratio=float(row["final_over_initial"]),
                late=float(row["late_window_l2_change"]),
                last=float(row["last_step_l2_change"]),
                z=float(row["final_peak_z_m"]),
                peak=float(row["final_peak"]),
            )
        )
    if active_rows:
        lines.append("")
        lines.append("## Active-Region Metrics")
        lines.append("")
        lines.append(
            "Temperatures and flows below use density weighting over cells above 0.1% of the same-snapshot density peak."
        )
        lines.append("")
        lines.append("| stage | sample | species | time [ms] | line density [m^-2] | peak n [m^-3] | active z [m] | <Tpar>_n [eV] | <Tper>_n [eV] | <u>_n [m/s] |")
        lines.append("|---|---|---|---:|---:|---:|---|---:|---:|---:|")
        for row in active_rows:
            lines.append(
                "| {stage} | {sample} | {species} | {time:.6f} | {line:.4e} | {peak:.4e} | {zmin:.3g}..{zmax:.3g} | {tpar:.4g} | {tper:.4g} | {u:.4g} |".format(
                    stage=row["stage"],
                    sample=row["sample"],
                    species=row["species"],
                    time=float(row["time_ms"]),
                    line=float(row["line_density_m-2"]),
                    peak=float(row["peak_density_m-3"]),
                    zmin=float(row["active_z_min_m"]),
                    zmax=float(row["active_z_max_m"]),
                    tpar=float(row.get("density_weighted_Tpar_m", math.nan)),
                    tper=float(row.get("density_weighted_Tper_m", math.nan)),
                    u=float(row.get("density_weighted_u_m", math.nan)),
                )
            )
    if particle_rows:
        lines.append("")
        lines.append("## Electron particle-energy metrics")
        lines.append("")
        lines.append("| sample | time [ms] | mean [eV] | p90 [eV] | p99 [eV] | p99.9 [eV] | frac >300 eV |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|")
        for row in particle_rows:
            lines.append(
                "| {sample} | {time:.6f} | {mean:.3g} | {p90:.3g} | {p99:.3g} | {p999:.3g} | {f300:.3e} |".format(
                    sample=row["snapshot_name"],
                    time=float(row["time_ms"]),
                    mean=float(row["mean_energy_eV"]),
                    p90=float(row["p90_energy_eV"]),
                    p99=float(row["p99_energy_eV"]),
                    p999=float(row["p999_energy_eV"]),
                    f300=float(row["frac_weight_E_gt_300eV"]),
                )
            )
    lines.append("")
    lines.append("## Figures")
    lines.append("")
    for fig in sorted(out_dir.glob("*.png")):
        lines.append(f"- `{fig.name}`")
    lines.append("")
    (out_dir / "mpex_ech0914_summary.md").write_text("\n".join(lines))


def main() -> int:
    args = parse_args()
    root = maybe_extract(args.path, args.extract_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    hdf5_dirs = discover_hdf5_dirs(root)
    stages = {key: load_stage(path) for key, path in hdf5_dirs.items()}

    moment_rows: list[dict[str, object]] = []
    active_rows: list[dict[str, object]] = []
    for stage in stages.values():
        moment_rows.extend(summary_rows(stage))
        active_rows.extend(active_region_rows(stage))
    write_csv(args.out_dir / "mpex_ech0914_moment_summary.csv", moment_rows)
    write_csv(args.out_dir / "mpex_ech0914_active_region_summary.csv", active_rows)

    plot_profiles(stages, args.out_dir, args.smooth_window)
    plot_density_zoom(stages, args.out_dir, args.smooth_window)
    for stage in stages.values():
        plot_time_maps(stage, args.out_dir, args.smooth_window, skip_initial=False)
        plot_time_maps(stage, args.out_dir, args.smooth_window, skip_initial=True)
    plot_fields(stages, args.out_dir, args.smooth_window)
    plot_eedf(stages, args.out_dir, args.electron_species, args.energy_max)
    if "ech" in stages:
        plot_energy_z(stages["ech"], args.out_dir, args.electron_species, args.energy_max)

    particle_rows: list[dict[str, object]] = []
    if "steady" in stages:
        particle_rows.extend(
            [
                particle_summary(stages["steady"], args.electron_species, "0", "steady start"),
                particle_summary(stages["steady"], args.electron_species, str(stages["steady"].times.size - 1), "steady final / ECH start"),
            ]
        )
    if "ech" in stages:
        particle_rows.append(
            particle_summary(stages["ech"], args.electron_species, str(stages["ech"].times.size - 1), "ECH final")
        )
    write_csv(args.out_dir / "mpex_ech0914_particle_energy_summary.csv", particle_rows)
    write_summary_markdown(args.out_dir, stages, moment_rows, active_rows, particle_rows)

    print(args.out_dir.resolve())
    print((args.out_dir / "mpex_ech0914_summary.md").resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
