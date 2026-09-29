#!/usr/bin/env python3
"""Compare the Nuclear Fusion Fig. 8 density with hybrid and kinetic inputs."""

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


def read_h5(file_path: Path, dataset: str) -> np.ndarray:
    executable = shutil.which("h5dump")
    if executable is None:
        raise RuntimeError("h5dump is required")
    text = subprocess.run(
        [executable, "-y", "-w", "0", "-d", dataset, str(file_path)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    match = re.search(r"DATA\s*\{(.*?)\}\s*\}\s*\}\s*$", text, re.S)
    if match is None:
        raise RuntimeError(f"Could not parse {dataset} from {file_path}")
    return np.fromstring(match.group(1).replace(",", " "), sep=" ")


def relative_l2(x: np.ndarray, values: np.ndarray, ref_x: np.ndarray, ref: np.ndarray) -> float:
    target = np.interp(x, ref_x, ref)
    return float(np.linalg.norm(values - target) / np.linalg.norm(target))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--kinetic-root", type=Path, default=None)
    parser.add_argument("--label", default="kinetic test")
    args = parser.parse_args()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    desktop = Path("/Users/78k/Desktop")
    reference_file = desktop / (
        "PICOS_paper_master_benchmark_20260924/"
        "PICOS_master_paper_fig8_9_profiles.npz"
    )
    hybrid_dir = desktop / "PICOS_fully_kinetic_matched_validation_20260925/hybrid"
    kinetic_root = args.kinetic_root.resolve() if args.kinetic_root else desktop / (
        "picos_kinetic_electron_eval/PICOS/validation/"
        "local_uncapped_steady_ech_validation_20260929_v7/profile_control_1us/"
        "picosFILES"
    )

    reference_data = np.load(reference_file)
    paper_x = reference_data["x_reference"]
    paper_density = reference_data["density_reference"]

    hybrid_x = read_h5(hybrid_dir / "main.h5", "/geometry/x_m")
    hybrid_profiles = np.asarray(
        [
            read_h5(hybrid_dir / "PARTICLES_FILE_0.h5", f"/{step}/ions/spp_1/n_m")
            for step in range(10, 23)
        ]
    )
    hybrid_density = np.mean(hybrid_profiles, axis=0)

    input_file = kinetic_root / "inputFiles/input_file.input"
    cv_ne = None
    for line in input_file.read_text().splitlines():
        if line.strip().startswith("CV_ne"):
            cv_ne = float(line.split()[1])
            break
    if cv_ne is None:
        raise RuntimeError("CV_ne was not found")
    profile_path = next((kinetic_root / "inputFiles").glob("*ne_norm_from_hybrid.txt"))
    kinetic_input_density = np.loadtxt(profile_path) * cv_ne
    kinetic_input_x = np.linspace(-2.0, 8.0, kinetic_input_density.size)

    kinetic_hdf = kinetic_root / "outputFiles/HDF5"
    kinetic_x = read_h5(kinetic_hdf / "main.h5", "/geometry/x_m")
    kinetic_initial_density = read_h5(
        kinetic_hdf / "PARTICLES_FILE_0.h5", "/0/ions/spp_2/n_m"
    )
    kinetic_final_density = read_h5(kinetic_hdf / "PARTICLES_FILE_0.h5", "/10/ions/spp_2/n_m")

    hybrid_error = relative_l2(hybrid_x, hybrid_density, paper_x, paper_density)
    input_error = relative_l2(
        kinetic_input_x, kinetic_input_density, paper_x, paper_density
    )
    final_error = relative_l2(
        kinetic_x, kinetic_final_density, paper_x, paper_density
    )
    initial_error = relative_l2(kinetic_x, kinetic_initial_density, paper_x, paper_density)

    fig, ax = plt.subplots(figsize=(11.5, 7.2))
    ax.plot(
        paper_x,
        paper_density / 1e20,
        "k--",
        lw=2.8,
        label="Nuclear Fusion Fig. 8(a), digitized",
    )
    ax.plot(
        hybrid_x,
        hybrid_density / 1e20,
        color="#228833",
        lw=2.5,
        label=f"validated hybrid steady average (relative L2={hybrid_error:.1%})",
    )
    ax.plot(
        kinetic_input_x,
        kinetic_input_density / 1e20,
        color="#4477aa",
        lw=2.2,
        label=f"profile supplied to {args.label} (relative L2={input_error:.1%})",
    )
    ax.plot(kinetic_x, kinetic_initial_density / 1e20, color="#aa4499", lw=1.5,
            label=f"{args.label} initial density (relative L2={initial_error:.1%})")
    ax.plot(
        kinetic_x,
        kinetic_final_density / 1e20,
        color="#cc3311",
        lw=2.0,
        label=f"{args.label} ECH-OFF final density (relative L2={final_error:.1%})",
    )
    ax.set_xlim(-2.0, 8.0)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel("z [m]")
    ax.set_ylabel("electron density [$10^{20}$ m$^{-3}$]")
    ax.set_title("Required density gate before kinetic ECH validation")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=10)
    fig.tight_layout()
    figure_path = out_dir / "nf_fig8_density_validation_gate.png"
    fig.savefig(figure_path, dpi=260)
    plt.close(fig)

    metrics_path = out_dir / "nf_fig8_density_validation_gate.txt"
    metrics_path.write_text(
        "\n".join(
            [
                f"hybrid_steady_average_relative_l2={hybrid_error:.9e}",
                f"hybrid_steady_average_peak_m3={np.max(hybrid_density):.9e}",
                f"kinetic_input_relative_l2={input_error:.9e}",
                f"kinetic_input_peak_m3={np.max(kinetic_input_density):.9e}",
                f"kinetic_initial_relative_l2={initial_error:.9e}",
                f"kinetic_initial_peak_m3={np.max(kinetic_initial_density):.9e}",
                f"ech_off_final_relative_l2={final_error:.9e}",
                f"ech_off_final_peak_m3={np.max(kinetic_final_density):.9e}",
            ]
        )
        + "\n"
    )
    print(figure_path)
    print(metrics_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
