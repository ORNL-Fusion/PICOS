#!/usr/bin/env python3
"""Plot matched PICOS++ ECH-OFF/ECH-ON profiles and electron EEDFs.

This reader intentionally uses the HDF5 command-line tools so it works on the
local validation machine without the Python h5py package.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


E_CHARGE = 1.602176634e-19
M_E = 9.1093837015e-31
COLORS = {"off": "#202020", "on": "#d95f02"}


def h5dump_values(file_path: Path, dataset: str) -> np.ndarray:
    executable = shutil.which("h5dump")
    if executable is None:
        raise RuntimeError("h5dump is required")
    result = subprocess.run(
        [executable, "-y", "-w", "0", "-d", dataset, str(file_path)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    match = re.search(r"DATA\s*\{(.*?)\}\s*\}\s*\}\s*$", result.stdout, re.S)
    if match is None:
        raise RuntimeError(f"Could not parse {dataset} from {file_path}")
    return np.fromstring(match.group(1).replace(",", " "), sep=" ", dtype=np.float64)


def smooth(values: np.ndarray, width: int = 3) -> np.ndarray:
    if width <= 1:
        return values
    kernel = np.ones(width, dtype=np.float64) / width
    return np.convolve(np.pad(values, width // 2, mode="edge"), kernel, mode="valid")


def final_step(hdf5_dir: Path) -> int:
    executable = shutil.which("h5ls")
    if executable is None:
        raise RuntimeError("h5ls is required")
    result = subprocess.run(
        [executable, str(hdf5_dir / "PARTICLES_FILE_0.h5")],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    steps = [int(line.split()[0]) for line in result.stdout.splitlines() if line.split()[0].isdigit()]
    return max(steps)


def read_profiles(hdf5_dir: Path, step: int) -> dict[str, np.ndarray]:
    main = hdf5_dir / "main.h5"
    particles = hdf5_dir / "PARTICLES_FILE_0.h5"
    fields = hdf5_dir / "FIELDS_FILE_0.h5"
    data = {"z": h5dump_values(main, "/geometry/x_m")}
    for species in ("spp_1", "spp_2"):
        for name in ("n_m", "Tpar_m", "Tper_m"):
            data[f"{species}_{name}"] = h5dump_values(
                particles, f"/{step}/ions/{species}/{name}"
            )
        data[f"{species}_u_m"] = h5dump_values(
            particles, f"/{step}/ions/{species}/u_m/x"
        )
    for name in ("BX_m", "EX_m"):
        data[name] = h5dump_values(fields, f"/{step}/fields/{name}/x")
    data["time"] = h5dump_values(fields, f"/{step}/time")
    return data


def read_electrons(hdf5_dir: Path, step: int) -> dict[str, np.ndarray]:
    chunks: dict[str, list[np.ndarray]] = {key: [] for key in ("z", "w", "vpar", "vper")}
    for particle_file in sorted(hdf5_dir.glob("PARTICLES_FILE_*.h5")):
        base = f"/{step}/ions/spp_2"
        z = h5dump_values(particle_file, f"{base}/X_p")
        w = h5dump_values(particle_file, f"{base}/a_p")
        velocity = h5dump_values(particle_file, f"{base}/V_p")
        if velocity.size != 2 * z.size:
            raise RuntimeError(f"Unexpected electron velocity shape in {particle_file}")
        velocity = velocity.reshape(2, z.size)
        valid = np.isfinite(z) & np.isfinite(w) & (w > 0.0)
        chunks["z"].append(z[valid])
        chunks["w"].append(w[valid])
        chunks["vpar"].append(velocity[0, valid])
        chunks["vper"].append(velocity[1, valid])
    result = {key: np.concatenate(values) for key, values in chunks.items()}
    result["eperp"] = 0.5 * M_E * result["vper"] ** 2 / E_CHARGE
    result["energy"] = (
        0.5 * M_E * (result["vpar"] ** 2 + result["vper"] ** 2) / E_CHARGE
    )
    return result


def eedf(energy: np.ndarray, weights: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    hist, _ = np.histogram(energy, bins=edges, weights=weights)
    widths = np.diff(edges)
    norm = np.sum(hist)
    density = hist / (norm * widths) if norm > 0.0 else np.zeros_like(hist)
    return 0.5 * (edges[:-1] + edges[1:]), density


def resonance_mean_profile(
    particles: dict[str, np.ndarray], z_edges: np.ndarray, quantity: str
) -> np.ndarray:
    numerator, _ = np.histogram(
        particles["z"], bins=z_edges, weights=particles["w"] * particles[quantity]
    )
    denominator, _ = np.histogram(particles["z"], bins=z_edges, weights=particles["w"])
    return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator > 0)


def style_axis(ax: plt.Axes) -> None:
    ax.grid(True, alpha=0.22)
    ax.tick_params(direction="in", top=True, right=True)


def plot_profiles(
    profiles: dict[str, dict[str, np.ndarray]], out_dir: Path, x1: float, x2: float
) -> Path:
    fig, axes = plt.subplots(3, 2, figsize=(14.5, 12.0), sharex=True)
    labels = {"off": "ECH OFF", "on": "ECH ON"}
    for key in ("off", "on"):
        p = profiles[key]
        z = p["z"]
        color = COLORS[key]
        axes[0, 0].plot(z, smooth(p["spp_2_n_m"]) / 1e20, color=color, lw=2.0, label=labels[key])
        axes[0, 0].plot(z, smooth(p["spp_1_n_m"]) / 1e20, color=color, lw=1.2, ls="--")
        axes[0, 1].plot(z, smooth(p["spp_2_Tpar_m"]), color=color, lw=2.0, label=labels[key])
        axes[1, 0].plot(z, smooth(p["spp_2_Tper_m"]), color=color, lw=2.0, label=labels[key])
        axes[1, 1].plot(z, smooth(p["spp_2_u_m"]) / 1e6, color=color, lw=2.0, label=labels[key])
        axes[2, 0].plot(z, smooth(p["EX_m"]), color=color, lw=1.8, label=labels[key])
        axes[2, 1].plot(z, smooth(p["spp_1_Tpar_m"]), color=color, lw=2.0, label=f"{labels[key]} $T_{{i\parallel}}$")
        axes[2, 1].plot(z, smooth(p["spp_1_Tper_m"]), color=color, lw=1.5, ls="--", label=f"{labels[key]} $T_{{i\perp}}$")

    axes[0, 0].set_title("Density (solid: electrons, dashed: D$^+$)")
    axes[0, 0].set_ylabel("$n$ [$10^{20}$ m$^{-3}$]")
    axes[0, 1].set_title("Electron parallel temperature")
    axes[0, 1].set_ylabel("$T_{e\parallel}$ [eV]")
    axes[1, 0].set_title("Electron perpendicular temperature")
    axes[1, 0].set_ylabel("$T_{e\perp}$ [eV]")
    axes[1, 1].set_title("Electron parallel flow")
    axes[1, 1].set_ylabel("$u_{e\parallel}$ [$10^6$ m/s]")
    axes[2, 0].set_title("Axial electric field")
    axes[2, 0].set_ylabel("$E_\parallel$ [V/m]")
    axes[2, 1].set_title("Ion temperatures")
    axes[2, 1].set_ylabel("$T_i$ [eV]")
    for ax in axes[-1, :]:
        ax.set_xlabel("z [m]")
    for ax in axes.ravel():
        ax.axvspan(x1, x2, color="#fdb863", alpha=0.18, lw=0)
        ax.set_xlim(-2.0, 8.0)
        style_axis(ax)
    for ax in axes.ravel():
        ax.legend(frameon=False, fontsize=9, ncol=2)
    fig.suptitle(
        "PICOS++ matched validation: ECH OFF vs ECH ON\n"
        "70 GHz, second harmonic, 300 kW; shaded region = 2.6–2.8 m",
        fontsize=16,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path = out_dir / "matched_ech_on_off_profiles.png"
    fig.savefig(path, dpi=240)
    plt.close(fig)
    return path


def plot_eedf_comparison(
    particles: dict[str, dict[str, np.ndarray]], out_dir: Path, x1: float, x2: float
) -> Path:
    e_edges = np.linspace(0.0, 300.0, 151)
    local_e_edges = np.linspace(0.0, 300.0, 61)
    map_e_edges = np.linspace(0.0, 300.0, 61)
    z_edges = np.linspace(-2.0, 8.0, 51)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    e_centers = 0.5 * (map_e_edges[:-1] + map_e_edges[1:])
    labels = {"off": "ECH OFF", "on": "ECH ON"}

    maps: dict[str, np.ndarray] = {}
    for key in ("off", "on"):
        p = particles[key]
        hist, _, _ = np.histogram2d(
            p["z"], p["energy"], bins=(z_edges, map_e_edges), weights=p["w"]
        )
        maps[key] = hist.T
    common_max = max(float(np.max(value)) for value in maps.values())

    fig = plt.figure(figsize=(15.5, 11.5))
    grid = fig.add_gridspec(2, 2, hspace=0.28, wspace=0.23)
    ax_global = fig.add_subplot(grid[0, 0])
    ax_local = fig.add_subplot(grid[0, 1])
    ax_off = fig.add_subplot(grid[1, 0])
    ax_on = fig.add_subplot(grid[1, 1], sharex=ax_off, sharey=ax_off)

    for key in ("off", "on"):
        p = particles[key]
        centers, density = eedf(p["energy"], p["w"], e_edges)
        ax_global.semilogy(centers, np.maximum(density, 1e-12), color=COLORS[key], lw=2.2, label=labels[key])
        local = (p["z"] >= x1) & (p["z"] <= x2)
        centers, density = eedf(p["energy"][local], p["w"][local], local_e_edges)
        ax_local.semilogy(centers, np.maximum(density, 1e-12), color=COLORS[key], lw=2.2, label=labels[key])

    for ax, title in ((ax_global, "Global electron energy distribution"), (ax_local, "EEDF in the ECH region (2.6–2.8 m)")):
        ax.set_xlabel("electron kinetic energy [eV]")
        ax.set_ylabel("weighted PDF [eV$^{-1}$]")
        ax.set_xlim(0.0, 300.0)
        ax.set_ylim(1e-6, 0.08)
        ax.set_title(title)
        ax.legend(frameon=False)
        style_axis(ax)

    image = None
    for key, ax in (("off", ax_off), ("on", ax_on)):
        normalized = np.log10(np.maximum(maps[key] / max(common_max, 1e-300), 1e-6))
        image = ax.pcolormesh(z_centers, e_centers, normalized, shading="auto", cmap="magma", vmin=-6.0, vmax=0.0)
        ax.axvspan(x1, x2, color="cyan", alpha=0.12, lw=0)
        ax.axvline(x1, color="cyan", ls="--", lw=1.2)
        ax.axvline(x2, color="cyan", ls="--", lw=1.2)
        ax.set_xlabel("z [m]")
        ax.set_ylabel("electron kinetic energy [eV]")
        ax.set_title(f"{labels[key]}: EEDF versus z")
        ax.set_xlim(-2.0, 8.0)
        ax.set_ylim(0.0, 300.0)
        style_axis(ax)
    if image is not None:
        cbar = fig.colorbar(image, ax=[ax_off, ax_on], location="bottom", shrink=0.78, pad=0.12)
        cbar.set_label("log$_{10}$ weighted counts, normalized to common maximum")

    fig.suptitle(
        "PICOS++ electron distribution validation: ECH OFF vs ECH ON\n"
        "Identical energy bins, weights, and color normalization",
        fontsize=16,
    )
    path = out_dir / "matched_ech_on_off_eedf.png"
    fig.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_local_energy_profiles(
    particles: dict[str, dict[str, np.ndarray]], out_dir: Path, x1: float, x2: float
) -> Path:
    z_edges = np.linspace(-2.0, 8.0, 101)
    z = 0.5 * (z_edges[:-1] + z_edges[1:])
    fig, axes = plt.subplots(2, 1, figsize=(12.5, 8.5), sharex=True)
    for key in ("off", "on"):
        color = COLORS[key]
        label = "ECH OFF" if key == "off" else "ECH ON"
        axes[0].plot(z, resonance_mean_profile(particles[key], z_edges, "eperp"), color=color, lw=2.1, label=label)
        axes[1].plot(z, resonance_mean_profile(particles[key], z_edges, "energy"), color=color, lw=2.1, label=label)
    axes[0].set_ylabel("weighted mean $E_\perp$ [eV]")
    axes[0].set_title("Particle-derived perpendicular electron energy")
    axes[1].set_ylabel("weighted mean total energy [eV]")
    axes[1].set_title("Particle-derived total electron energy")
    axes[1].set_xlabel("z [m]")
    for ax in axes:
        ax.axvspan(x1, x2, color="#fdb863", alpha=0.2, lw=0)
        ax.set_xlim(-2.0, 8.0)
        ax.legend(frameon=False)
        style_axis(ax)
    fig.suptitle("Spatial localization of the ECH response", fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path = out_dir / "matched_ech_on_off_particle_energy_profiles.png"
    fig.savefig(path, dpi=240)
    plt.close(fig)
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--x1", type=float, default=2.6)
    parser.add_argument("--x2", type=float, default=2.8)
    args = parser.parse_args()
    root = args.root.resolve()
    out_dir = (args.out_dir or root / "plots").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    dirs = {
        "off": root / "profile_control_1us/picosFILES/outputFiles/HDF5",
        "on": root / "profile_ech_1us/picosFILES/outputFiles/HDF5",
    }
    steps = {key: final_step(value) for key, value in dirs.items()}
    profiles = {key: read_profiles(dirs[key], steps[key]) for key in dirs}
    particles = {key: read_electrons(dirs[key], steps[key]) for key in dirs}
    paths = [
        plot_profiles(profiles, out_dir, args.x1, args.x2),
        plot_eedf_comparison(particles, out_dir, args.x1, args.x2),
        plot_local_energy_profiles(particles, out_dir, args.x1, args.x2),
    ]
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
