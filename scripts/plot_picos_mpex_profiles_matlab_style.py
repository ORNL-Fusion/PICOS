#!/usr/bin/env python3
"""Plot PICOS++ MPEX profiles in the ProtoLite MATLAB-preview style.

This is a Python/HDF5 equivalent of the data path used by
GetDataHDF5 1.m and PICOS_MPEX_PreviewData 1.m:

* read main.h5 for geometry and MPI layout,
* read FIELDS_FILE_0.h5 for mesh fields,
* read PARTICLES_FILE_0.h5 for reduced mesh moments,
* read all PARTICLES_FILE_*.h5 only for particle-state diagnostics.

The original MATLAB preview assumes one ion species and accidentally reads
spp_1 mesh moments for every species.  This script keeps that data layout but
uses the selected species name, which matters for kinetic-electron runs.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
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
AMU = 1.66053906660e-27
M_E = 9.1093837015e-31

SPECIES_LABELS = {
    "spp_1": "D+",
    "spp_2": "electron",
}


@dataclass
class PicosCase:
    name: str
    hdf5_dir: Path
    input_dir: Path | None
    input_params: dict[str, str]
    ion_params: dict[str, str]
    z_m: np.ndarray
    times_s: np.ndarray
    fields: dict[str, np.ndarray]
    moments: dict[str, dict[str, np.ndarray]]
    species: list[str]


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


def parse_key_value_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path or not path.exists():
        return values
    for raw in path.read_text(errors="ignore").splitlines():
        line = raw.split("//", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 2:
            values[parts[0]] = parts[1]
    return values


def find_case_input_dir(hdf5_dir: Path) -> Path | None:
    for parent in [hdf5_dir, *hdf5_dir.parents]:
        candidate = parent / "picosFILES" / "inputFiles"
        if candidate.is_dir():
            return candidate
        if parent.name == "outputFiles":
            candidate = parent.parent / "inputFiles"
            if candidate.is_dir():
                return candidate
    return None


def find_input_file(input_dir: Path | None, case_name: str) -> Path | None:
    if input_dir is None:
        return None
    exact = input_dir / f"input_file_{case_name}.input"
    if exact.exists():
        return exact
    matches = sorted(input_dir.glob("input_file_*.input"))
    return matches[0] if matches else None


def find_ion_file(input_dir: Path | None, case_name: str) -> Path | None:
    if input_dir is None:
        return None
    exact = input_dir / f"ions_properties_{case_name}.ion"
    if exact.exists():
        return exact
    matches = sorted(input_dir.glob("ions_properties*.ion"))
    return matches[0] if matches else None


def load_profile(input_dir: Path | None, file_name: str | None) -> np.ndarray | None:
    if input_dir is None or not file_name or file_name == "none":
        return None
    path = input_dir / file_name
    if not path.exists():
        return None
    values = np.loadtxt(path, dtype=np.float64)
    return np.asarray(values, dtype=np.float64).reshape(-1)


def profile_axis(params: dict[str, str], n: int) -> np.ndarray:
    xmin = float(params.get("LX_min", "0.0"))
    xmax = float(params.get("LX_max", str(n)))
    return np.linspace(xmin, xmax, n, dtype=np.float64)


def find_hdf5_dirs(path: Path, extract_dir: Path | None) -> list[Path]:
    path = path.expanduser()
    if path.is_file() and path.suffixes[-2:] == [".tar", ".gz"]:
        if extract_dir is None:
            extract_dir = Path("/private/tmp") / f"picos_mpex_profiles_{path.stem.replace('.tar', '')}"
        extract_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(path, "r:gz") as tar:
            tar.extractall(extract_dir)
        return find_hdf5_dirs(extract_dir, None)
    if path.is_file() and path.name == "main.h5":
        return [path.parent.resolve()]
    if path.is_dir() and (path / "main.h5").exists():
        return [path.resolve()]
    if path.is_dir():
        return sorted({candidate.parent.resolve() for candidate in path.rglob("HDF5/main.h5")})
    raise FileNotFoundError(path)


def load_case(hdf5_dir: Path) -> PicosCase:
    case_name = hdf5_dir.parent.name
    input_dir = find_case_input_dir(hdf5_dir)
    input_file = find_input_file(input_dir, case_name)
    ion_file = find_ion_file(input_dir, case_name)
    input_params = parse_key_value_file(input_file) if input_file else {}
    ion_params = parse_key_value_file(ion_file) if ion_file else {}

    with h5py.File(hdf5_dir / "main.h5", "r") as handle:
        z_m = read_1d(handle, "geometry/x_m")
        if z_m is None:
            z_m = read_1d(handle, "geometry/xAxis")
        if z_m is None:
            raise KeyError(f"{hdf5_dir / 'main.h5'} has no geometry/x_m or geometry/xAxis")

    fields: dict[str, np.ndarray] = {}
    times_s: np.ndarray | None = None
    with h5py.File(hdf5_dir / "FIELDS_FILE_0.h5", "r") as handle:
        steps = numeric_keys(handle)
        if not steps:
            raise ValueError(f"No numeric snapshots in {hdf5_dir / 'FIELDS_FILE_0.h5'}")
        times_s = np.asarray([read_scalar(handle, f"{step}/time") for step in steps], dtype=np.float64)
        for name in ("BX_m", "dBX_m", "ddBX_m", "EX_m", "Phi_m"):
            cols = []
            for step in steps:
                values = read_1d(handle, f"{step}/fields/{name}/x")
                if values is None:
                    cols = []
                    break
                cols.append(values)
            if cols:
                fields[name] = np.column_stack(cols)

    p0 = hdf5_dir / "PARTICLES_FILE_0.h5"
    if not p0.exists():
        raise FileNotFoundError(p0)

    with h5py.File(p0, "r") as handle:
        psteps = numeric_keys(handle)
        if psteps != steps:
            steps = psteps
            times_s = np.asarray([read_scalar(handle, f"{step}/time") for step in steps], dtype=np.float64)
        species = sorted(handle[f"{steps[0]}/ions"].keys(), key=lambda item: int(item.split("_")[-1]))
        moments: dict[str, dict[str, np.ndarray]] = {}
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

    n = min([z_m.size] + [arr.shape[0] for group in moments.values() for arr in group.values()] + [arr.shape[0] for arr in fields.values()])
    return PicosCase(
        name=case_name,
        hdf5_dir=hdf5_dir,
        input_dir=input_dir,
        input_params=input_params,
        ion_params=ion_params,
        z_m=z_m[:n],
        times_s=np.asarray(times_s, dtype=np.float64),
        fields={key: values[:n, :] for key, values in fields.items()},
        moments={spp: {key: values[:n, :] for key, values in group.items()} for spp, group in moments.items()},
        species=species,
    )


def moving_average(values: np.ndarray, window: int, axis: int = 0) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if window <= 1:
        return arr
    kernel = np.ones(window, dtype=np.float64) / float(window)
    return np.apply_along_axis(lambda x: np.convolve(x, kernel, mode="same"), axis, arr)


def finite_window_mean(values: np.ndarray, indices: Iterable[int]) -> np.ndarray:
    idx = list(indices)
    if not idx:
        idx = [values.shape[1] - 1]
    subset = np.asarray(values[:, idx], dtype=np.float64)
    valid = np.isfinite(subset)
    counts = np.count_nonzero(valid, axis=1)
    sums = np.nansum(np.where(valid, subset, 0.0), axis=1)
    out = np.full(subset.shape[0], np.nan, dtype=np.float64)
    np.divide(sums, counts, out=out, where=counts > 0)
    return out


def selected_time_indices(n_times: int) -> list[int]:
    if n_times <= 1:
        return [0]
    return sorted(set([0, max(0, n_times // 2), n_times - 1]))


def mask_temperature(temp: np.ndarray, density: np.ndarray, floor_fraction: float) -> np.ndarray:
    out = np.asarray(temp, dtype=np.float64).copy()
    peak = np.nanmax(density) if np.any(np.isfinite(density)) else 0.0
    if peak > 0.0:
        out[density < peak * floor_fraction] = np.nan
    return out


def interp_profile_to_mesh(profile: np.ndarray | None, z_mesh: np.ndarray, params: dict[str, str]) -> np.ndarray | None:
    if profile is None:
        return None
    zp = profile_axis(params, profile.size)
    return np.interp(z_mesh, zp, profile)


def scaled_overlay(ax: plt.Axes, z: np.ndarray, values: np.ndarray | None, target_max: float, label: str, color: str, style: str = "-") -> None:
    if values is None or not np.any(np.isfinite(values)):
        return
    vmax = float(np.nanmax(np.abs(values)))
    if vmax <= 0.0 or not np.isfinite(vmax) or target_max <= 0.0 or not np.isfinite(target_max):
        return
    ax.plot(z, values * (target_max / vmax), style, color=color, lw=1.8, alpha=0.85, label=label)


def plot_source_and_inputs(case: PicosCase, out_dir: Path) -> None:
    params = case.input_params
    ion = case.ion_params
    z = case.z_m
    ne = interp_profile_to_mesh(load_profile(case.input_dir, ion.get("IC_densityFraction_fileName_1")), z, params)
    te = interp_profile_to_mesh(load_profile(case.input_dir, params.get("IC_Te_fileName")), z, params)
    ti = interp_profile_to_mesh(load_profile(case.input_dir, ion.get("IC_Tpar_fileName_1")), z, params)
    source = interp_profile_to_mesh(load_profile(case.input_dir, params.get("pairSource_fileName") or ion.get("BC_G_fileName_1")), z, params)
    b = interp_profile_to_mesh(load_profile(case.input_dir, params.get("IC_BX_fileName")), z, params)

    fig, axes = plt.subplots(2, 1, figsize=(12.5, 8.6), sharex=True)
    ax = axes[0]
    if ne is not None:
        ax.plot(z, ne, "k", lw=2.4, label="IC density fraction")
    if source is not None:
        ax.plot(z, source, color="#d62728", lw=2.2, label="pair source shape")
    if te is not None:
        ax.plot(z, te, color="#1f77b4", lw=2.0, label="Te fraction")
    if ti is not None:
        ax.plot(z, ti, color="#9467bd", lw=1.8, label="Ti fraction")
    ax.set_ylabel("normalized profile")
    ax.set_title("Input density, temperature, and source profiles")
    ax.grid(True, alpha=0.24)
    ax.legend(frameon=False, ncol=2)

    ax = axes[1]
    if "BX_m" in case.fields:
        ax.plot(z, case.fields["BX_m"][:, 0], "k", lw=2.4, label="B from output")
    if b is not None:
        b0 = float(params.get("IC_BX", "1.0"))
        ax.plot(z, b0 * b, color="#d62728", lw=1.8, ls="--", label="B input profile")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("B [T]")
    ax.set_xlim(float(np.nanmin(z)), float(np.nanmax(z)))
    ax.grid(True, alpha=0.24)
    ax.legend(frameon=False)

    fig.tight_layout()
    fig.savefig(out_dir / f"{case.name}_input_profiles.png", dpi=220)
    plt.close(fig)


def plot_profile_dashboard(case: PicosCase, out_dir: Path, smooth_window: int, final_window: int, density_floor: float) -> None:
    z = case.z_m
    idxs = selected_time_indices(case.times_s.size)
    final_start = max(0, case.times_s.size - final_window)
    final_idx = range(final_start, case.times_s.size)

    fig, axes = plt.subplots(3, 2, figsize=(16.0, 13.0), sharex=True)
    axes = axes.reshape(3, 2)
    colors = ["0.35", "#1f77b4", "#d62728"]
    labels = [f"{case.times_s[i] * 1.0e3:.3g} ms" for i in idxs]

    for col, spp in enumerate(case.species[:2]):
        label = SPECIES_LABELS.get(spp, spp)
        ax = axes[0, col]
        density = case.moments[spp].get("n_m")
        if density is not None:
            density_curves = []
            for color, tlabel, ii in zip(colors, labels, idxs):
                curve = moving_average(density[:, ii], smooth_window)
                density_curves.append((ii, curve))
                ax.plot(z, curve, color=color, lw=2.0, label=tlabel)
            final_n = moving_average(finite_window_mean(density, final_idx), smooth_window)
            density_curves.append((-1, final_n))
            ax.plot(z, final_n, color="k", lw=2.8, ls="--", label="final avg")
            scale_curves = [curve for ii, curve in density_curves if ii != 0]
            scale_values = np.concatenate([curve[np.isfinite(curve)] for curve in scale_curves if np.any(np.isfinite(curve))])
            ymax = float(np.nanpercentile(scale_values, 99.5)) * 1.35 if scale_values.size else 1.0
            if ymax > 0.0 and np.isfinite(ymax):
                ax.set_ylim(bottom=0.0, top=ymax)
                if np.any(np.isfinite(density_curves[0][1])) and np.nanmax(density_curves[0][1]) > 2.0 * ymax:
                    ax.text(0.02, 0.88, "t=0 off scale", transform=ax.transAxes, fontsize=9, color="0.35")
            if col == 0:
                source = interp_profile_to_mesh(
                    load_profile(case.input_dir, case.input_params.get("pairSource_fileName") or case.ion_params.get("BC_G_fileName_1")),
                    z,
                    case.input_params,
                )
                scaled_overlay(ax, z, source, ax.get_ylim()[1] * 0.95, "source shape", "#2ca02c", ":")
        ax.set_title(f"{label} density, linear scale")
        ax.set_ylabel("n [m^-3]")
        ax.grid(True, alpha=0.24)
        ax.legend(frameon=False, fontsize=9)

        ax = axes[1, col]
        tpar = case.moments[spp].get("Tpar_m")
        tper = case.moments[spp].get("Tper_m")
        density = case.moments[spp].get("n_m")
        if tpar is not None and density is not None:
            nbar = finite_window_mean(density, final_idx)
            tpar_bar = mask_temperature(finite_window_mean(tpar, final_idx), nbar, density_floor)
            ax.plot(z, moving_average(tpar_bar, smooth_window), "k", lw=2.2, label="Tpar")
        if tper is not None and density is not None:
            nbar = finite_window_mean(density, final_idx)
            tper_bar = mask_temperature(finite_window_mean(tper, final_idx), nbar, density_floor)
            ax.plot(z, moving_average(tper_bar, smooth_window), color="#d62728", lw=2.0, label="Tper")
        ax.set_title(f"{label} final-window temperature")
        ax.set_ylabel("T [eV]")
        ax.grid(True, alpha=0.24)
        ax.legend(frameon=False)

        ax = axes[2, col]
        u = case.moments[spp].get("u_m")
        if u is not None:
            ubar = moving_average(finite_window_mean(u, final_idx), smooth_window)
            if spp == "spp_1" and "spp_2" in case.moments:
                te = case.moments["spp_2"].get("Tpar_m")
                ti = case.moments["spp_1"].get("Tpar_m")
                den = case.moments["spp_1"].get("n_m")
                if te is not None and ti is not None and den is not None:
                    te_bar = mask_temperature(finite_window_mean(te, final_idx), finite_window_mean(den, final_idx), density_floor)
                    ti_bar = mask_temperature(finite_window_mean(ti, final_idx), finite_window_mean(den, final_idx), density_floor)
                    m_ion = 2.0 * AMU
                    cs = np.sqrt(E_CHARGE * np.maximum(te_bar + 3.0 * ti_bar, 0.0) / m_ion)
                    mach = ubar / np.maximum(cs, 1.0e-300)
                    ax.plot(z, mach, "k", lw=2.2, label="D+ Mach")
                    ax.set_ylabel("U / Cs")
                else:
                    ax.plot(z, ubar, "k", lw=2.2, label="u")
                    ax.set_ylabel("u [m/s]")
            else:
                ax.plot(z, ubar, "k", lw=2.2, label="u")
                ax.set_ylabel("u [m/s]")
        ax.set_title(f"{label} final-window flow")
        ax.set_xlabel("z [m]")
        ax.grid(True, alpha=0.24)
        ax.legend(frameon=False)

    for ax in axes.reshape(-1):
        ax.set_xlim(float(np.nanmin(z)), float(np.nanmax(z)))
    fig.suptitle(f"{case.name}: ProtoLite-style PICOS++ profile preview")
    fig.tight_layout()
    fig.savefig(out_dir / f"{case.name}_profile_dashboard.png", dpi=220)
    plt.close(fig)


def plot_time_maps(case: PicosCase, out_dir: Path, smooth_window: int, skip_initial: bool) -> None:
    suffix = "no_t0" if skip_initial and case.times_s.size > 1 else "with_t0"
    tidx = np.arange(case.times_s.size)
    if skip_initial and case.times_s.size > 1:
        tidx = tidx[1:]
    z = case.z_m
    t_ms = case.times_s[tidx] * 1.0e3

    panels: list[tuple[str, str, str, np.ndarray]] = []
    for spp in case.species[:2]:
        label = SPECIES_LABELS.get(spp, spp)
        density = case.moments[spp].get("n_m")
        tpar = case.moments[spp].get("Tpar_m")
        if density is not None:
            panels.append((f"{label} density", "n [m^-3]", "viridis", moving_average(density[:, tidx], smooth_window)))
        if tpar is not None and density is not None:
            temp = tpar[:, tidx].copy()
            for jj in range(temp.shape[1]):
                temp[:, jj] = mask_temperature(temp[:, jj], density[:, tidx[jj]], 1.0e-4)
            panels.append((f"{label} Tpar", "Tpar [eV]", "magma", moving_average(temp, smooth_window)))
    if not panels:
        return

    rows = int(math.ceil(len(panels) / 2))
    fig, axes = plt.subplots(rows, 2, figsize=(15.5, 4.9 * rows), sharex=True, squeeze=False)
    for ax, (title, cbar_label, cmap, values) in zip(axes.reshape(-1), panels):
        finite = values[np.isfinite(values)]
        vmax = float(np.nanpercentile(finite, 99.5)) if finite.size else 1.0
        if vmax <= 0.0 or not np.isfinite(vmax):
            vmax = 1.0
        im = ax.pcolormesh(z, t_ms, values.T, shading="auto", cmap=cmap, vmin=0.0, vmax=vmax)
        ax.set_title(f"{title} ({'excluding t=0' if suffix == 'no_t0' else 'including t=0'})")
        ax.set_ylabel("time [ms]")
        ax.grid(True, alpha=0.18)
        cbar = fig.colorbar(im, ax=ax, pad=0.01)
        cbar.set_label(cbar_label)
    for ax in axes.reshape(-1)[len(panels):]:
        ax.axis("off")
    for ax in axes[-1, :]:
        ax.set_xlabel("z [m]")
    fig.tight_layout()
    fig.savefig(out_dir / f"{case.name}_time_maps_{suffix}.png", dpi=220)
    plt.close(fig)


def plot_matlab_mesh_surfaces(case: PicosCase, out_dir: Path, smooth_window: int) -> None:
    species = [spp for spp in case.species[:2] if "n_m" in case.moments.get(spp, {})]
    if not species:
        return
    fig = plt.figure(figsize=(15.5, 6.2 * len(species)))
    t_grid, z_grid = np.meshgrid(case.times_s * 1.0e3, case.z_m)
    for ii, spp in enumerate(species, start=1):
        values = moving_average(case.moments[spp]["n_m"], smooth_window)
        ax = fig.add_subplot(len(species), 1, ii, projection="3d")
        ax.plot_surface(t_grid, z_grid, values, cmap="viridis", linewidth=0, antialiased=False)
        ax.set_title(f"{SPECIES_LABELS.get(spp, spp)} density mesh(t,z,n), MATLAB preview style")
        ax.set_xlabel("time [ms]", labelpad=8)
        ax.set_ylabel("z [m]", labelpad=8)
        ax.set_zlabel("n [m^-3]", labelpad=8)
        ax.view_init(elev=24.0, azim=-135.0)
    fig.tight_layout()
    fig.savefig(out_dir / f"{case.name}_density_mesh_surfaces.png", dpi=220)
    plt.close(fig)


def plot_field_dashboard(case: PicosCase, out_dir: Path, smooth_window: int) -> None:
    z = case.z_m
    final_start = max(0, case.times_s.size - 5)
    final_idx = range(final_start, case.times_s.size)
    fig, axes = plt.subplots(2, 2, figsize=(15.0, 9.0), sharex=True)
    axes = axes.reshape(-1)

    if "BX_m" in case.fields:
        b = finite_window_mean(case.fields["BX_m"], final_idx)
        axes[0].plot(z, moving_average(b, smooth_window), "k", lw=2.3)
        axes[0].set_ylabel("B [T]")
        axes[0].set_title("final-window B")
    if "EX_m" in case.fields:
        e = finite_window_mean(case.fields["EX_m"], final_idx)
        axes[1].plot(z, moving_average(e, smooth_window), "k", lw=2.3)
        axes[1].set_ylabel("E_parallel [V/m]")
        axes[1].set_title("final-window E")
    if "Phi_m" in case.fields:
        phi = finite_window_mean(case.fields["Phi_m"], final_idx)
        axes[2].plot(z, moving_average(phi, smooth_window), "k", lw=2.3)
        axes[2].set_ylabel("Phi [V]")
        axes[2].set_title("final-window potential")
    if "BX_m" in case.fields and "EX_m" in case.fields:
        e = moving_average(case.fields["EX_m"], smooth_window)
        finite = e[np.isfinite(e)]
        vmax = float(np.nanpercentile(np.abs(finite), 99.0)) if finite.size else 1.0
        vmax = max(vmax, 1.0)
        im = axes[3].pcolormesh(z, case.times_s * 1.0e3, e.T, shading="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
        axes[3].set_title("E_parallel time map")
        axes[3].set_ylabel("time [ms]")
        cbar = fig.colorbar(im, ax=axes[3], pad=0.01)
        cbar.set_label("E_parallel [V/m]")
    for ax in axes:
        ax.grid(True, alpha=0.24)
        ax.set_xlim(float(np.nanmin(z)), float(np.nanmax(z)))
    for ax in axes[2:]:
        ax.set_xlabel("z [m]")
    fig.suptitle(f"{case.name}: field profiles")
    fig.tight_layout()
    fig.savefig(out_dir / f"{case.name}_field_dashboard.png", dpi=220)
    plt.close(fig)


def line_integral(z: np.ndarray, values: np.ndarray) -> float:
    mask = np.isfinite(z) & np.isfinite(values)
    if np.count_nonzero(mask) < 2:
        return math.nan
    return float(np.trapz(values[mask], z[mask]))


def summarize(case: PicosCase, out_dir: Path, final_window: int) -> None:
    rows: list[dict[str, object]] = []
    final_start = max(0, case.times_s.size - final_window)
    previous_start = max(0, final_start - final_window)
    for spp in case.species:
        label = SPECIES_LABELS.get(spp, spp)
        density = case.moments[spp].get("n_m")
        if density is None:
            continue
        initial = density[:, 0]
        final = finite_window_mean(density, range(final_start, case.times_s.size))
        previous = finite_window_mean(density, range(previous_start, final_start)) if final_start > previous_start else initial
        peak_idx = int(np.nanargmax(final)) if np.any(np.isfinite(final)) else 0
        rows.append(
            {
                "case": case.name,
                "species": label,
                "quantity": "density",
                "n_snapshots": case.times_s.size,
                "final_time_ms": float(case.times_s[-1] * 1.0e3),
                "initial_line_integral_m-2": line_integral(case.z_m, initial),
                "final_window_line_integral_m-2": line_integral(case.z_m, final),
                "final_to_initial_line_integral_ratio": line_integral(case.z_m, final) / max(line_integral(case.z_m, initial), 1.0e-300),
                "previous_to_final_l2_change": float(np.linalg.norm(final - previous) / max(np.linalg.norm(previous), 1.0e-300)),
                "last_step_l2_change": float(np.linalg.norm(density[:, -1] - density[:, -2]) / max(np.linalg.norm(density[:, -2]), 1.0e-300)) if density.shape[1] > 1 else math.nan,
                "final_window_peak_m-3": float(np.nanmax(final)),
                "final_window_peak_z_m": float(case.z_m[peak_idx]),
            }
        )

    input_dir = case.input_dir
    source = load_profile(input_dir, case.input_params.get("pairSource_fileName") or case.ion_params.get("BC_G_fileName_1"))
    ne = load_profile(input_dir, case.ion_params.get("IC_densityFraction_fileName_1"))
    with (out_dir / f"{case.name}_profile_summary.csv").open("w", newline="") as fp:
        keys = list(rows[0].keys()) if rows else ["case"]
        writer = csv.DictWriter(fp, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    lines = [f"# PICOS++ MPEX Profile Preview: {case.name}", ""]
    lines.append(f"HDF5: `{case.hdf5_dir}`")
    if input_dir:
        lines.append(f"Input files: `{input_dir}`")
    lines.append(f"Snapshots: {case.times_s.size}, final saved time: {case.times_s[-1] * 1.0e3:.6g} ms")
    if source is not None:
        zsrc = profile_axis(case.input_params, source.size)
        lines.append(f"Pair source peak: z = {zsrc[int(np.nanargmax(source))]:.4g} m")
    if ne is not None:
        zne = profile_axis(case.input_params, ne.size)
        lines.append(f"IC density-fraction peak: z = {zne[int(np.nanargmax(ne))]:.4g} m")
    lines.extend(["", "| species | final/inital line density | late-window L2 change | last-step L2 change | final peak z [m] | final peak n [m^-3] |", "|---|---:|---:|---:|---:|---:|"])
    for row in rows:
        lines.append(
            "| {species} | {ratio:.4g} | {window:.4g} | {last:.4g} | {z:.4g} | {peak:.4e} |".format(
                species=row["species"],
                ratio=float(row["final_to_initial_line_integral_ratio"]),
                window=float(row["previous_to_final_l2_change"]),
                last=float(row["last_step_l2_change"]),
                z=float(row["final_window_peak_z_m"]),
                peak=float(row["final_window_peak_m-3"]),
            )
        )
    lines.append("")
    lines.append("Interpretation note: temperature/flow profiles are masked in cells where the same-species density is below the configured density floor, because moment temperatures in nearly empty cells are not physically meaningful.")
    (out_dir / f"{case.name}_profile_summary.md").write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path, help="HDF5 dir, output root, or .tar.gz archive.")
    parser.add_argument("--out-dir", type=Path, default=Path("validation/mpex_profiles_matlab_style"))
    parser.add_argument("--extract-dir", type=Path, default=None)
    parser.add_argument("--case-filter", default="", help="Only analyze HDF5 paths containing this substring.")
    parser.add_argument("--smooth-window", type=int, default=10)
    parser.add_argument("--final-window", type=int, default=5)
    parser.add_argument("--temperature-density-floor", type=float, default=1.0e-3)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    hdf5_dirs: list[Path] = []
    for path in args.paths:
        hdf5_dirs.extend(find_hdf5_dirs(path, args.extract_dir))
    hdf5_dirs = sorted(set(hdf5_dirs))
    if args.case_filter:
        hdf5_dirs = [path for path in hdf5_dirs if args.case_filter in str(path)]
    if not hdf5_dirs:
        print("No HDF5/main.h5 directories found.")
        return 2

    seen_names: set[str] = set()
    for hdf5_dir in hdf5_dirs:
        case = load_case(hdf5_dir)
        if case.name in seen_names:
            print(f"Skipping duplicate case name: {case.name}")
            continue
        seen_names.add(case.name)
        case_out = args.out_dir / case.name
        case_out.mkdir(parents=True, exist_ok=True)
        plot_source_and_inputs(case, case_out)
        plot_profile_dashboard(case, case_out, args.smooth_window, args.final_window, args.temperature_density_floor)
        plot_time_maps(case, case_out, args.smooth_window, skip_initial=False)
        plot_time_maps(case, case_out, args.smooth_window, skip_initial=True)
        plot_matlab_mesh_surfaces(case, case_out, args.smooth_window)
        plot_field_dashboard(case, case_out, args.smooth_window)
        summarize(case, case_out, args.final_window)
        print(case_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
