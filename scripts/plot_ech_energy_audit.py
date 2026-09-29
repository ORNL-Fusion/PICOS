#!/usr/bin/env python3
"""Audit ECH energy conservation and ON/OFF electron heating histories."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_matched_ech_validation import final_step, h5dump_values

E_CHARGE = 1.602176634e-19


def history(hdf: Path) -> dict[str, np.ndarray]:
    main = hdf / "main.h5"
    ncp = float(h5dump_values(main, "/ions/spp_2/NCP")[0])
    mass = float(h5dump_values(main, "/ions/spp_2/M")[0])
    result = {key: [] for key in ("time", "energy", "mean", "region", "eperp_region")}
    for step in range(final_step(hdf) + 1):
        total_energy = total_number = region_energy = region_perp = region_number = 0.0
        for particle_file in sorted(hdf.glob("PARTICLES_FILE_*.h5")):
            base = f"/{step}/ions/spp_2"
            x = h5dump_values(particle_file, f"{base}/X_p")
            w = h5dump_values(particle_file, f"{base}/a_p")
            velocity = h5dump_values(particle_file, f"{base}/V_p").reshape(2, x.size)
            eperp = 0.5*mass*velocity[1]**2
            energy = 0.5*mass*np.sum(velocity**2, axis=0)
            real_weight = ncp*w
            selected = (x >= 2.6) & (x <= 2.8)
            total_energy += np.sum(real_weight*energy)
            total_number += np.sum(real_weight)
            region_energy += np.sum(real_weight[selected]*energy[selected])
            region_perp += np.sum(real_weight[selected]*eperp[selected])
            region_number += np.sum(real_weight[selected])
        result["time"].append(float(h5dump_values(hdf / "PARTICLES_FILE_0.h5", f"/{step}/time")[0]))
        result["energy"].append(total_energy)
        result["mean"].append(total_energy/total_number/E_CHARGE)
        result["region"].append(region_energy/region_number/E_CHARGE)
        result["eperp_region"].append(region_perp/region_number/E_CHARGE)
    return {key: np.asarray(value) for key, value in result.items()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    out = root / "plots"
    out.mkdir(exist_ok=True)
    paths = {
        "off": root / "profile_control_1us/picosFILES/outputFiles/HDF5",
        "on": root / "profile_ech_1us/picosFILES/outputFiles/HDF5",
    }
    data = {key: history(path) for key, path in paths.items()}
    boundary = {}
    for key, path in paths.items():
        particle_file = path / "PARTICLES_FILE_0.h5"
        boundary[key] = {
            name: np.asarray([
                float(h5dump_values(particle_file, f"/{step}/boundary/{name}")[0])
                for step in range(final_step(path) + 1)
            ])
            for name in ("cumulativeElectronEnergyLost", "cumulativeElectronEnergyInjected")
        }
    on_file = paths["on"] / "PARTICLES_FILE_0.h5"
    steps = range(final_step(paths["on"]) + 1)
    absorbed = np.asarray([
        float(h5dump_values(on_file, f"/{step}/rf/electron/cumulativeAbsorbedEnergy")[0])
        for step in steps
    ])
    target = np.asarray([
        float(h5dump_values(on_file, f"/{step}/rf/electron/cumulativeTargetEnergy")[0])
        for step in steps
    ])

    time_us = data["on"]["time"]*1e6
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    axes[0, 0].plot(time_us, data["off"]["energy"], label="ECH OFF", color="k")
    axes[0, 0].plot(time_us, data["on"]["energy"], label="ECH ON", color="#d95f02")
    axes[0, 0].set_ylabel("represented electron energy [J]")
    axes[0, 1].plot(time_us, data["on"]["energy"] - data["off"]["energy"], color="#0072b2", label="ON − OFF")
    extra_net_boundary_loss = (
        boundary["on"]["cumulativeElectronEnergyLost"]
        - boundary["on"]["cumulativeElectronEnergyInjected"]
        - boundary["off"]["cumulativeElectronEnergyLost"]
        + boundary["off"]["cumulativeElectronEnergyInjected"]
    )
    stored_plus_boundary = data["on"]["energy"] - data["off"]["energy"] + extra_net_boundary_loss
    axes[0, 1].plot(time_us, stored_plus_boundary, color="#009e73",
                    label="ON − OFF + extra net boundary loss")
    axes[0, 1].plot(time_us, absorbed, "--", color="#d95f02", label="RF absorbed")
    axes[0, 1].plot(time_us, target, ":", color="k", label="RF target")
    axes[0, 1].set_ylabel("energy difference [J]")
    for key, label, color in (("off", "ECH OFF", "k"), ("on", "ECH ON", "#d95f02")):
        axes[1, 0].plot(time_us, data[key]["mean"], color=color, label=label)
        axes[1, 1].plot(time_us, data[key]["eperp_region"], color=color, label=label)
    axes[1, 0].set_ylabel("global mean electron energy [eV]")
    axes[1, 1].set_ylabel(r"mean $E_\perp$ in 2.6–2.8 m [eV]")
    for ax in axes.ravel():
        ax.grid(alpha=0.25)
        ax.legend(frameon=False)
    axes[1, 0].set_xlabel("time [µs]")
    axes[1, 1].set_xlabel("time [µs]")
    fig.suptitle("PICOS++ ECH energy audit: 300 kW second-harmonic electron heating")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    figure = out / "ech_energy_audit.png"
    fig.savefig(figure, dpi=240)
    plt.close(fig)

    final = -1
    metrics = out / "ech_energy_audit.txt"
    metrics.write_text(
        f"duration_s={data['on']['time'][final]:.9e}\n"
        f"rf_target_J={target[final]:.9e}\n"
        f"rf_absorbed_J={absorbed[final]:.9e}\n"
        f"rf_absorbed_to_target={absorbed[final]/target[final]:.9e}\n"
        f"particle_energy_on_minus_off_J={data['on']['energy'][final]-data['off']['energy'][final]:.9e}\n"
        f"global_mean_off_eV={data['off']['mean'][final]:.9e}\n"
        f"global_mean_on_eV={data['on']['mean'][final]:.9e}\n"
        f"region_Eperp_off_eV={data['off']['eperp_region'][final]:.9e}\n"
        f"region_Eperp_on_eV={data['on']['eperp_region'][final]:.9e}\n"
        f"electron_boundary_loss_off_J={boundary['off']['cumulativeElectronEnergyLost'][final]:.9e}\n"
        f"electron_boundary_loss_on_J={boundary['on']['cumulativeElectronEnergyLost'][final]:.9e}\n"
        f"electron_injected_energy_off_J={boundary['off']['cumulativeElectronEnergyInjected'][final]:.9e}\n"
        f"electron_injected_energy_on_J={boundary['on']['cumulativeElectronEnergyInjected'][final]:.9e}\n"
        f"extra_net_boundary_loss_on_minus_off_J={extra_net_boundary_loss[final]:.9e}\n"
        f"stored_plus_extra_net_boundary_loss_J={stored_plus_boundary[final]:.9e}\n"
    )
    print(figure)
    print(metrics)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
