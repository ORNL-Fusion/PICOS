#!/usr/bin/env python3
"""Analyze PICOS ECH/ICH HDF5 outputs.

This replaces the hard-coded ProtoLite post-processing scripts with a command
line workflow. It reads PICOS HDF5 output, summarizes ion/electron heating,
plots field and plasma profiles, computes force-balance terms, and optionally
builds velocity-space histograms with the correct guiding-center
``2*pi*v_perp`` Jacobian.
"""

from __future__ import annotations

import argparse
import csv
import math
import warnings
from pathlib import Path
from typing import Any

import numpy as np

from picos_hdf5 import (
    C_LIGHT,
    E_CHARGE,
    PicosRun,
    central_diff,
    kinetic_energy_eV,
    read_field_series,
    read_species_mesh_series,
    read_species_particle_state,
    rolling_mean,
    species_metadata,
    velocity_components,
    weighted_mean,
    weighted_percentile,
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _as_float(value: Any, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _rf_block(main: dict[str, Any], z: float) -> tuple[str, dict[str, Any], int]:
    rf = main.get("rf", {})
    if z < 0.0:
        return "ECH", rf.get("electron", {}), int(rf.get("heatElectrons", 0))
    return "ICH", rf.get("ion", {}), int(rf.get("heatIons", 0))


def _run_label(run: PicosRun) -> str:
    if run.hdf5_dir.name == "HDF5":
        return run.hdf5_dir.parent.name
    return run.hdf5_dir.name


def _select_time_window(n_steps: int, last: int) -> slice:
    count = max(1, min(last, n_steps))
    return slice(n_steps - count, n_steps)


def _field_value(fields: dict[str, np.ndarray], name: str, default_shape: tuple[int, int]) -> np.ndarray:
    if name in fields:
        return fields[name]
    return np.full(default_shape, np.nan)


def _nanmean(values: np.ndarray, axis: int | None = None) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(values, axis=axis)


def _nanmin(values: np.ndarray) -> float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return float(np.nanmin(values))


def _nanmax(values: np.ndarray) -> float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return float(np.nanmax(values))


def _nanrms(values: np.ndarray) -> float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return float(np.sqrt(np.nanmean(values * values)))


def _mesh_mask(x_m: np.ndarray, x_window: tuple[float, float] | None) -> np.ndarray:
    if x_window is None:
        return np.isfinite(x_m)
    x1, x2 = x_window
    lo = min(x1, x2)
    hi = max(x1, x2)
    return np.isfinite(x_m) & (x_m >= lo) & (x_m <= hi)


def _particle_mask(state: dict[str, np.ndarray], x_window: tuple[float, float] | None) -> np.ndarray:
    if "V_p" not in state:
        return np.array([], dtype=bool)
    n_particles = velocity_components(state["V_p"]).shape[1]
    mask = np.ones(n_particles, dtype=bool)
    if x_window is not None and "X_p" in state:
        x1, x2 = x_window
        lo = min(x1, x2)
        hi = max(x1, x2)
        mask &= (state["X_p"] >= lo) & (state["X_p"] <= hi)
    return mask


def particle_temperature_summary(
    state: dict[str, np.ndarray], mass_kg: float, x_window: tuple[float, float] | None
) -> dict[str, float]:
    if "V_p" not in state:
        return {}

    velocity = velocity_components(state["V_p"])
    mask = _particle_mask(state, x_window)
    if not np.any(mask):
        return {
            "particle_count_selected": 0.0,
            "particle_Tpar_eV": math.nan,
            "particle_Tper_eV": math.nan,
            "particle_upar_m_s": math.nan,
        }

    velocity = velocity[:, mask]
    weights = state.get("a_p")
    if weights is not None:
        weights = np.asarray(weights, dtype=float)[mask]

    vpar = velocity[0]
    vperp = np.hypot(velocity[1], velocity[2]) if velocity.shape[0] > 2 else np.abs(velocity[1])
    upar = weighted_mean(vpar, weights)
    if weights is None:
        varpar = float(np.mean((vpar - upar) ** 2))
        vperp2 = float(np.mean(vperp * vperp))
    else:
        valid = np.isfinite(weights) & (weights > 0.0)
        varpar = float(np.average((vpar[valid] - upar) ** 2, weights=weights[valid])) if np.any(valid) else math.nan
        vperp2 = float(np.average(vperp[valid] * vperp[valid], weights=weights[valid])) if np.any(valid) else math.nan

    return {
        "particle_count_selected": float(vpar.size),
        "particle_Tpar_eV": mass_kg * varpar / E_CHARGE,
        "particle_Tper_eV": mass_kg * vperp2 / (2.0 * E_CHARGE),
        "particle_upar_m_s": upar,
    }


def summarize_species(
    run: PicosRun,
    species: str,
    times: np.ndarray,
    fields: dict[str, np.ndarray],
    last: int,
    x_window: tuple[float, float] | None,
    particle_stride: int,
    max_particles: int | None,
) -> dict[str, Any]:
    meta = species_metadata(run.main, species)
    z = meta["Z"]
    channel, rf, heat_flag = _rf_block(run.main, z)
    time_slice = _select_time_window(times.size, last)
    x_m = run.x_m
    mask_x = _mesh_mask(x_m, x_window)

    mesh = read_species_mesh_series(run.hdf5_dir, species, ("n_m", "Tpar_m", "Tper_m", "u_m"), steps=run.steps)
    row: dict[str, Any] = {
        "run": _run_label(run),
        "hdf5_dir": str(run.hdf5_dir),
        "species": species,
        "channel": channel,
        "Z": z,
        "mass_kg": meta["M"],
        "charge_C": meta["Q"],
        "rf_enabled_for_channel": heat_flag,
        "rf_power_W": _as_float(rf.get("Prf")),
        "rf_frequency_Hz": _as_float(rf.get("freq")),
        "rf_harmonic": int(_as_float(rf.get("n_harmonic"), 0.0)),
        "rf_x1_m": _as_float(rf.get("x1")),
        "rf_x2_m": _as_float(rf.get("x2")),
        "rf_efield_mode": int(_as_float(rf.get("eFieldMode"), 0.0)),
        "rf_efield_amplitude_Vm": _as_float(rf.get("eFieldAmplitude")),
        "fieldSolveModel": int(_as_float(run.main.get("fieldSolveModel"), 0.0)),
        "relativisticElectrons": int(_as_float(run.main.get("relativisticElectrons"), 0.0)),
        "t_start_s": float(times[0]) if times.size else math.nan,
        "t_end_s": float(times[-1]) if times.size else math.nan,
        "n_outputs": int(times.size),
    }

    for key, values in mesh.items():
        selected = values[mask_x, time_slice]
        row[f"{key}_mean_last"] = float(_nanmean(selected)) if selected.size else math.nan
        row[f"{key}_max_last"] = _nanmax(selected) if selected.size else math.nan
        row[f"{key}_min_last"] = _nanmin(selected) if selected.size else math.nan

    shape = (x_m.size, times.size)
    ex = _field_value(fields, "EX_m", shape)
    bx = _field_value(fields, "BX_m", shape)
    row["EX_m_rms_last_Vm"] = _nanrms(ex[mask_x, time_slice])
    row["BX_m_mean_last_T"] = float(_nanmean(bx[mask_x, time_slice]))

    final_step = run.steps[-1]
    state = read_species_particle_state(
        run.hdf5_dir,
        species,
        final_step,
        stride=max(1, particle_stride),
        max_particles=max_particles,
    )
    if "V_p" in state:
        rel = row["relativisticElectrons"] == 1 and z < 0.0
        energy = kinetic_energy_eV(state["V_p"], meta["M"], relativistic=rel)
        weights = state.get("a_p")
        row["particle_energy_mean_eV"] = weighted_mean(energy, weights)
        row["particle_energy_p50_eV"] = weighted_percentile(energy, 50, weights)
        row["particle_energy_p95_eV"] = weighted_percentile(energy, 95, weights)
        row["particle_energy_p99_eV"] = weighted_percentile(energy, 99, weights)
        row["particle_energy_max_eV"] = float(np.nanmax(energy)) if energy.size else math.nan
        row["particle_max_speed_over_c"] = float(
            np.nanmax(np.sqrt(np.sum(velocity_components(state["V_p"]) ** 2, axis=0))) / C_LIGHT
        )
        row.update(particle_temperature_summary(state, meta["M"], x_window))

    return row


def write_summary_csv(rows: list[dict[str, Any]], out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "ech_ich_summary.csv"
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def plot_profiles(
    run: PicosRun,
    times: np.ndarray,
    fields: dict[str, np.ndarray],
    species_list: list[str],
    last: int,
    smooth: int,
    out_dir: Path,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    x_m = run.x_m
    time_slice = _select_time_window(times.size, last)
    n_cols = 2
    fig, axes = plt.subplots(3, n_cols, figsize=(12, 11), sharex=True)
    axes = axes.ravel()

    if "BX_m" in fields:
        axes[0].plot(x_m, np.nanmean(fields["BX_m"][:, time_slice], axis=1), "k", lw=1.8, label="B")
    axes[0].set_ylabel("B [T]")
    axes[0].set_title("Magnetic field")
    axes[0].grid(True, alpha=0.25)

    if "EX_m" in fields:
        ex = rolling_mean(np.nanmean(fields["EX_m"][:, time_slice], axis=1), smooth, axis=0)
        axes[1].plot(x_m, ex, "tab:red", lw=1.8, label="E_x")
    axes[1].set_ylabel("E_x [V/m]")
    axes[1].set_title("Electric field")
    axes[1].grid(True, alpha=0.25)

    for species in species_list:
        mesh = read_species_mesh_series(run.hdf5_dir, species, ("n_m", "Tpar_m", "Tper_m", "u_m"), steps=run.steps)
        meta = species_metadata(run.main, species)
        channel = "e-" if meta["Z"] < 0.0 else "ion"
        if "n_m" in mesh:
            n_prof = rolling_mean(_nanmean(mesh["n_m"][:, time_slice], axis=1), smooth, axis=0)
            axes[2].plot(x_m, n_prof, lw=1.5, label=f"{species} {channel}")
        if "Tpar_m" in mesh:
            tpar = rolling_mean(_nanmean(mesh["Tpar_m"][:, time_slice], axis=1), smooth, axis=0)
            axes[3].plot(x_m, tpar, lw=1.5, label=f"{species} Tpar")
        if "Tper_m" in mesh:
            tper = rolling_mean(_nanmean(mesh["Tper_m"][:, time_slice], axis=1), smooth, axis=0)
            axes[3].plot(x_m, tper, lw=1.5, ls="--", label=f"{species} Tper")
        if "u_m" in mesh:
            u_prof = rolling_mean(_nanmean(mesh["u_m"][:, time_slice], axis=1), smooth, axis=0)
            axes[4].plot(x_m, u_prof, lw=1.5, label=species)

    axes[2].set_ylabel("n [m^-3]")
    axes[2].set_title("Density")
    axes[2].grid(True, alpha=0.25)
    axes[2].legend(fontsize=8)

    axes[3].set_ylabel("T [eV]")
    axes[3].set_title("Parallel/perpendicular temperature")
    axes[3].grid(True, alpha=0.25)
    axes[3].legend(fontsize=8)

    axes[4].set_ylabel("u_parallel [m/s]")
    axes[4].set_title("Parallel flow")
    axes[4].grid(True, alpha=0.25)
    axes[4].legend(fontsize=8)

    if "Phi_m" in fields:
        phi = rolling_mean(_nanmean(fields["Phi_m"][:, time_slice], axis=1), smooth, axis=0)
        axes[5].plot(x_m, phi, "tab:purple", lw=1.8)
    axes[5].set_ylabel("phi [V]")
    axes[5].set_title("Electrostatic potential")
    axes[5].grid(True, alpha=0.25)

    for ax in axes[-2:]:
        ax.set_xlabel("x [m]")

    fig.suptitle(f"PICOS ECH/ICH profiles: {_run_label(run)}", fontsize=14)
    fig.tight_layout()
    path = out_dir / "ech_ich_profiles.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_time_traces(
    run: PicosRun,
    times: np.ndarray,
    species_list: list[str],
    x_window: tuple[float, float] | None,
    out_dir: Path,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x_m = run.x_m
    mask_x = _mesh_mask(x_m, x_window)
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)

    for species in species_list:
        mesh = read_species_mesh_series(run.hdf5_dir, species, ("n_m", "Tpar_m", "Tper_m"), steps=run.steps)
        if "n_m" in mesh:
            axes[0].plot(times, _nanmean(mesh["n_m"][mask_x, :], axis=0), lw=1.5, label=species)
        if "Tpar_m" in mesh:
            axes[1].plot(times, _nanmean(mesh["Tpar_m"][mask_x, :], axis=0), lw=1.5, label=f"{species} Tpar")
        if "Tper_m" in mesh:
            axes[2].plot(times, _nanmean(mesh["Tper_m"][mask_x, :], axis=0), lw=1.5, label=f"{species} Tper")

    axes[0].set_ylabel("mean n [m^-3]")
    axes[1].set_ylabel("mean Tpar [eV]")
    axes[2].set_ylabel("mean Tper [eV]")
    axes[2].set_xlabel("time [s]")
    for ax in axes:
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle(f"PICOS ECH/ICH time traces: {_run_label(run)}", fontsize=14)
    fig.tight_layout()
    path = out_dir / "ech_ich_time_traces.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_force_balance(
    run: PicosRun,
    times: np.ndarray,
    fields: dict[str, np.ndarray],
    species: str,
    last: int,
    smooth: int,
    out_dir: Path,
) -> Path | None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mesh = read_species_mesh_series(run.hdf5_dir, species, ("n_m", "Tpar_m", "Tper_m", "u_m"), steps=run.steps)
    if not all(key in mesh for key in ("n_m", "Tpar_m", "Tper_m", "u_m")) or "BX_m" not in fields:
        return None

    meta = species_metadata(run.main, species)
    x_m = run.x_m
    dx = float(np.mean(np.diff(x_m)))
    time_slice = _select_time_window(times.size, last)
    n = mesh["n_m"]
    tpar = mesh["Tpar_m"]
    tper = mesh["Tper_m"]
    u = mesh["u_m"]
    b = fields["BX_m"]

    p_par = E_CHARGE * tpar * n
    p_per = E_CHARGE * tper * n
    p_ke = meta["M"] * n * u * u
    d_p_par_dx = central_diff(p_par, dx, axis=0)
    d_b_dx = central_diff(b, dx, axis=0)
    f_ke_x = b * central_diff(p_ke / np.maximum(b, 1.0e-30), dx, axis=0)
    f_par = d_p_par_dx
    f_mag = (p_per - p_par) / np.maximum(b, 1.0e-30) * d_b_dx

    f_ke_x = rolling_mean(_nanmean(f_ke_x[:, time_slice], axis=1), smooth, axis=0)
    f_par = rolling_mean(_nanmean(f_par[:, time_slice], axis=1), smooth, axis=0)
    f_mag = rolling_mean(_nanmean(f_mag[:, time_slice], axis=1), smooth, axis=0)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x_m, f_ke_x, "k", lw=1.8, label="kinetic-gradient")
    ax.plot(x_m, f_par, "tab:blue", lw=1.8, label="parallel pressure")
    ax.plot(x_m, f_mag, "tab:red", lw=1.8, label="mirror")
    ax.axhline(0.0, color="0.35", lw=0.8)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("force density [N m^-3]")
    ax.set_title(f"Force balance terms: {_run_label(run)} {species}")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    path = out_dir / f"force_balance_{species}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_velocity_space(
    run: PicosRun,
    species: str,
    step: str,
    x_window: tuple[float, float] | None,
    bins: int,
    particle_stride: int,
    max_particles: int | None,
    out_dir: Path,
) -> Path | None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    state = read_species_particle_state(
        run.hdf5_dir,
        species,
        step,
        stride=max(1, particle_stride),
        max_particles=max_particles,
    )
    if "V_p" not in state:
        return None

    velocity = velocity_components(state["V_p"])
    mask = _particle_mask(state, x_window)
    if not np.any(mask):
        return None
    velocity = velocity[:, mask]
    weights = state.get("a_p")
    if weights is None:
        weights = np.ones(velocity.shape[1])
    else:
        weights = np.asarray(weights, dtype=float)[mask]

    vpar = velocity[0]
    vperp = np.hypot(velocity[1], velocity[2]) if velocity.shape[0] > 2 else np.abs(velocity[1])
    finite = np.isfinite(vpar) & np.isfinite(vperp) & np.isfinite(weights) & (weights > 0.0)
    if np.count_nonzero(finite) < 10:
        return None
    vpar = vpar[finite]
    vperp = vperp[finite]
    weights = weights[finite]

    vpar_limit = np.nanpercentile(np.abs(vpar), 99.5)
    vperp_limit = np.nanpercentile(vperp, 99.5)
    if not np.isfinite(vpar_limit) or vpar_limit <= 0.0:
        vpar_limit = np.nanmax(np.abs(vpar))
    if not np.isfinite(vperp_limit) or vperp_limit <= 0.0:
        vperp_limit = np.nanmax(vperp)
    vpar_edges = np.linspace(-1.05 * vpar_limit, 1.05 * vpar_limit, bins + 1)
    vperp_edges = np.linspace(0.0, 1.05 * vperp_limit, bins + 1)

    counts, _, _ = np.histogram2d(vpar, vperp, bins=(vpar_edges, vperp_edges), weights=weights)
    dvpar = np.diff(vpar_edges)[:, None]
    dvperp = np.diff(vperp_edges)[None, :]
    total_weight = np.sum(weights)
    reduced_pdf = counts / max(total_weight, 1.0e-300) / np.maximum(dvpar * dvperp, 1.0e-300)

    vperp_centers = 0.5 * (vperp_edges[:-1] + vperp_edges[1:])
    jacobian = 2.0 * np.pi * np.maximum(vperp_centers[None, :], 0.5 * np.diff(vperp_edges)[0])
    full_gyro_pdf = reduced_pdf / jacobian

    meta = species_metadata(run.main, species)
    t_ref = max(1.0e-12, _as_float(run.main.get("Te"), 1.0))
    if meta["Z"] > 0.0:
        block = run.main.get("ions", {}).get(species, {})
        t_ref = max(1.0e-12, _as_float(block.get("Tpar"), t_ref))
    v_ref = math.sqrt(2.0 * E_CHARGE * t_ref / meta["M"]) if meta["M"] > 0.0 else 1.0

    x_centers = 0.5 * (vpar_edges[:-1] + vpar_edges[1:]) / v_ref
    y_centers = vperp_centers / v_ref
    floor = 1.0e-300

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
    im0 = axes[0].pcolormesh(
        x_centers,
        y_centers,
        np.log10(np.maximum(reduced_pdf.T, floor)),
        shading="auto",
        cmap="magma",
    )
    axes[0].set_title("reduced g(vpar,vperp)")
    axes[0].set_xlabel("v_parallel / v_ref")
    axes[0].set_ylabel("v_perp / v_ref")
    fig.colorbar(im0, ax=axes[0], label="log10 g")

    im1 = axes[1].pcolormesh(
        x_centers,
        y_centers,
        np.log10(np.maximum(full_gyro_pdf.T, floor)),
        shading="auto",
        cmap="viridis",
    )
    axes[1].set_title("full f = g/(2*pi*vperp)")
    axes[1].set_xlabel("v_parallel / v_ref")
    fig.colorbar(im1, ax=axes[1], label="log10 f")

    for ax in axes:
        ax.grid(False)
        ax.set_aspect("auto")

    fig.suptitle(f"Velocity-space diagnostic: {_run_label(run)} {species} step {step}", fontsize=13)
    fig.tight_layout()
    path = out_dir / f"velocity_space_{species}_step{step}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def write_notes(out_dir: Path, rows: list[dict[str, Any]], artifacts: list[Path]) -> Path:
    path = out_dir / "ech_ich_analysis.md"
    lines = [
        "# PICOS ECH/ICH Analysis",
        "",
        "This analysis reads PICOS HDF5 output directly and uses the species charge sign to label ECH and ICH channels.",
        "",
        "Moment convention:",
        "",
        "- PICOS mesh moments are direct particle-weight sums over the sampled guiding-center ensemble.",
        "- `P_perp = m n <v_perp^2>/2`; density and pressure deposition are not divided by `v_perp`.",
        "- The `1/v_perp` factor appears only when reconstructing the full gyrotropic distribution from a `(v_parallel,v_perp)` histogram: `f = g/(2*pi*v_perp)`.",
        "",
        "Summary rows:",
        "",
    ]
    for row in rows:
        lines.append(
            f"- `{row['run']}` `{row['species']}` {row['channel']}: "
            f"RF={row.get('rf_enabled_for_channel')} Prf={row.get('rf_power_W')} W, "
            f"Tpar_mean={row.get('Tpar_m_mean_last', math.nan):.6g} eV, "
            f"Tper_mean={row.get('Tper_m_mean_last', math.nan):.6g} eV"
        )
    lines.extend(["", "Artifacts:", ""])
    for artifact in artifacts:
        lines.append(f"- `{artifact.name}`")
    path.write_text("\n".join(lines) + "\n")
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", help="PICOS run tag, run directory, or direct HDF5 directory")
    parser.add_argument("--out-dir", type=Path, default=None, help="Directory for CSV/PNG/Markdown artifacts")
    parser.add_argument("--species", nargs="*", default=None, help="Species names such as spp_1 spp_2; default: all")
    parser.add_argument("--last", type=int, default=5, help="Average this many final output frames")
    parser.add_argument("--smooth", type=int, default=7, help="Spatial smoothing window for line plots")
    parser.add_argument("--x-window", nargs=2, type=float, default=None, metavar=("X1", "X2"))
    parser.add_argument("--particle-stride", type=int, default=1, help="Read every Nth output particle for diagnostics")
    parser.add_argument("--max-particles", type=int, default=None, help="Limit particles per species after striding")
    parser.add_argument("--velocity-hist", action="store_true", help="Write reduced and gyrotropic velocity histograms")
    parser.add_argument("--hist-bins", type=int, default=120)
    args = parser.parse_args()

    repo_root = _repo_root()
    run = PicosRun.open(args.output, repo_root=repo_root)
    out_dir = args.out_dir
    if out_dir is None:
        out_dir = repo_root / "validation" / "ech_ich_analysis" / _run_label(run)
    out_dir.mkdir(parents=True, exist_ok=True)

    species_list = args.species if args.species else run.species
    x_window = tuple(args.x_window) if args.x_window is not None else None
    times, fields = read_field_series(run.hdf5_dir, variables=("BX_m", "EX_m", "Phi_m", "dBX_m", "ddBX_m"), steps=run.steps)

    rows = [
        summarize_species(
            run,
            species,
            times,
            fields,
            last=args.last,
            x_window=x_window,
            particle_stride=args.particle_stride,
            max_particles=args.max_particles,
        )
        for species in species_list
    ]

    artifacts: list[Path] = []
    artifacts.append(write_summary_csv(rows, out_dir))
    artifacts.append(plot_profiles(run, times, fields, species_list, args.last, args.smooth, out_dir))
    artifacts.append(plot_time_traces(run, times, species_list, x_window, out_dir))
    for species in species_list:
        force_path = plot_force_balance(run, times, fields, species, args.last, args.smooth, out_dir)
        if force_path is not None:
            artifacts.append(force_path)
        if args.velocity_hist:
            velocity_path = plot_velocity_space(
                run,
                species,
                run.steps[-1],
                x_window,
                args.hist_bins,
                args.particle_stride,
                args.max_particles,
                out_dir,
            )
            if velocity_path is not None:
                artifacts.append(velocity_path)
    artifacts.append(write_notes(out_dir, rows, artifacts))

    print(f"Read: {run.hdf5_dir}")
    for artifact in artifacts:
        print(f"Wrote: {artifact}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
