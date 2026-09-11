#!/usr/bin/env python3
"""Validate PICOS++ source build-up before RF/ECH studies.

The validation uses fixed computational markers with zero initial physical
weight.  In explicit pair-source mode those zero-weight markers are activated
from a warm source centered at z=0 using the requested pair source rate.  The
expected qualitative behavior is:

* n_i and n_e are initially zero in every cell.
* density appears first near the source and then fills the open domain.
* finite-temperature cells stay close to the source birth temperature.
* the pair-source path injects equal-and-opposite charge pairs at the requested
  rate while recording particle, energy, and momentum source terms.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import textwrap
import time
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PICOS_ROOT = Path(__file__).resolve().parents[1]
SPECIES_LABEL = {"spp_1": "D+", "spp_2": "e-"}


@dataclass(frozen=True)
class SourceCase:
    tag: str
    pair_source: bool


def write_vector(path: Path, values: np.ndarray) -> None:
    np.savetxt(path, np.asarray(values, dtype=float), fmt="%.16e")


def input_deck(tag: str, pair_source: bool, source_rate: float, sim_gyroperiods: float, output_gyroperiods: float) -> str:
    pair_switch = 1 if pair_source else 0
    return textwrap.dedent(
        f"""
        // PICOS++ particle-source build-up validation.
        // Purpose: validate source-driven density/temperature build-up before ECH.
        // The computational markers are initialized with zero physical weight
        // through IC_weightScale=0.  The warm z=0 source then fills the domain.
        // =============================================================================
        mpisForFields               1
        quietStart                  1
        IC_velocityDistributionModel 0
        IC_randomSeed               424242
        IC_weightScale              0.0
        numberOfParticleSpecies     2
        numberOfTracerSpecies       0
        advanceParticleMethod       1

        // Characteristic values:
        // =============================================================================
        CV_ne                       1.0000000000000000e+19
        CV_Te                       5.0000000000000000e+03
        CV_B                        1.0000000000000001e-01
        CV_Tpar                     5.0000000000000000e+03
        CV_Tper                     5.0000000000000000e+03

        // Simulation time is in reference D+ gyroperiod units.
        // =============================================================================
        DTc                         1.0000000000000001e-01
        simulationTime              {sim_gyroperiods:.16e}

        // Switches:
        // =============================================================================
        SW_EfieldSolve              0
        SW_fieldSolveModel          0
        SW_electronGyroTimeStepLimiter 0
        SW_electronPlasmaTimeStepLimiter 0
        SW_BfieldSolve              0
        SW_Collisions               0
        CollOperType                2
        SW_collisionConservationProjection 0
        collisionRandomSeed         314159
        SW_RFheating                0
        SW_RFheatingIons            0
        SW_RFheatingElectrons       0
        SW_relativisticElectrons    0
        SW_pairSource               {pair_switch}
        SW_advancePos               1
        SW_linearSolve              0

        // Restart controls:
        // =============================================================================
        restart_enabled             0
        restart_path                none
        restart_snapshot            -1
        restart_continueTime        0
        restart_particleFilePrefix  PARTICLES_FILE_
        restart_fieldsFilePrefix    FIELDS_FILE_

        // Magnetic field initial conditions:
        // =============================================================================
        IC_uniformBfield            1
        IC_BX                       1.0000000000000001e-01
        IC_BY                       0.0
        IC_BZ                       0.0
        IC_BX_NX                    80
        IC_BX_fileName              {tag}_one.txt
        IC_phiLeft                  0.0
        IC_phiRight                 0.0
        Poisson_BCModel             1
        Poisson_sheathCoefficient   3.0
        ReformulatedPoisson_lambda  -1.0
        ReformulatedPoisson_quasiNeutral 0

        // Geometry:
        // =============================================================================
        dp                          5.0000000000000003e-02
        r1                          0.0
        r2                          2.0000000000000000e-02
        LX_min                      -2.0000000000000001e-01
        LX_max                      2.0000000000000001e-01

        // Electron fluid/profile initial conditions:
        // This profile is not used as a source of physical particles here; kinetic
        // electrons are species 2 and also start with IC_weightScale=0.
        // =============================================================================
        IC_ne                       1.0000000000000000e+19
        IC_Te                       5.0000000000000000e+03
        IC_Te_NX                    80
        IC_Te_fileName              {tag}_one.txt

        // Coupled electron-ion source:
        // pairSource_weightMode=1 uses pairSource_rate to set injected weights.
        // =============================================================================
        pairSource_ionSpecies       1
        pairSource_electronSpecies  2
        pairSource_rate             {source_rate:.16e}
        pairSource_mean_x           0.0
        pairSource_sigma_x          4.0000000000000001e-02
        pairSource_Ti_birth         5.0000000000000000e+03
        pairSource_Te_birth         5.0000000000000000e+03
        pairSource_Ei_birth         0.0
        pairSource_Ee_birth         0.0
        pairSource_eta_i            7.8539816339744828e-01
        pairSource_eta_e            7.8539816339744828e-01
        pairSource_positionMode     0
        pairSource_weightMode       1
        pairSource_fileName         {tag}_one.txt
        pairSource_NS               80
        pairSource_maxParticleWeight 1000

        // RF operator, disabled for source validation:
        // =============================================================================
        RF_ion_Prf                      0.0
        RF_ion_n_harmonic               1
        RF_ion_freq                     1.0
        RF_ion_x1                       -2.0000000000000001e-01
        RF_ion_x2                       2.0000000000000001e-01
        RF_ion_t_ON                     0.0
        RF_ion_t_OFF                    {sim_gyroperiods:.16e}
        RF_ion_kpar                     0.0
        RF_ion_kper                     0.0
        RF_ion_handedness               -1
        RF_ion_EfieldMode               1
        RF_ion_resonanceMode            0
        RF_ion_EfieldAmplitude          0.0
        RF_ion_maxEnergyGainFraction    0.0
        RF_ion_maxParticleEnergy        0.0
        RF_ion_maxVelocityFractionC     0.0
        RF_ion_Prf_fileName             {tag}_one.txt
        RF_ion_Prf_NS                   80

        RF_electron_Prf                      0.0
        RF_electron_n_harmonic               1
        RF_electron_freq                     1.0
        RF_electron_x1                       -2.0000000000000001e-01
        RF_electron_x2                       2.0000000000000001e-01
        RF_electron_t_ON                     0.0
        RF_electron_t_OFF                    {sim_gyroperiods:.16e}
        RF_electron_kpar                     0.0
        RF_electron_kper                     0.0
        RF_electron_handedness               -1
        RF_electron_EfieldMode               1
        RF_electron_resonanceMode            0
        RF_electron_EfieldAmplitude          0.0
        RF_electron_maxEnergyGainFraction    0.0
        RF_electron_maxParticleEnergy        0.0
        RF_electron_maxVelocityFractionC     0.0
        RF_electron_Prf_fileName             {tag}_one.txt
        RF_electron_Prf_NS                   80

        // Output variables:
        // =============================================================================
        outputCadence               {output_gyroperiods:.16e}
        outputs_variables           {{X_p,V_p,a_p,BX_p,n_m,Tpar_m,Tper_m,u_m}}

        // Data smoothing:
        // =============================================================================
        smoothingParameter          0.0
        filtersPerIterationFields   0
        filtersPerIterationIons     0
        """
    ).strip() + "\n"


def ion_deck(tag: str, npc: int, source_rate: float, source_sigma: float) -> str:
    return textwrap.dedent(
        f"""
        // PICOS++ source build-up validation species deck.
        // Both kinetic species start with zero physical marker weight.  In
        // pair-source mode they are fueled from z=0 using pairSource_rate.
        // In legacy mode each species uses its own BC_G_* boundary source.
        // =============================================================================
        SPECIES1                      1
        NPC1                          {npc}
        pctSupPartOutput1             100.0
        Z1                            1
        M1                            2.0000000000000000e+00

        IC_type_1                     1
        IC_Tper_1                     5.0000000000000000e+03
        IC_Tper_fileName_1            {tag}_one.txt
        IC_Tper_NX_1                  80
        IC_Tpar_1                     5.0000000000000000e+03
        IC_Tpar_fileName_1            {tag}_one.txt
        IC_Tpar_NX_1                  80
        IC_densityFraction_1          1.0
        IC_densityFraction_fileName_1 {tag}_one.txt
        IC_densityFraction_NX_1       80
        IC_weightScale_1              0.0

        BC_type_1                     1
        BC_T_1                        5.0000000000000000e+03
        BC_E_1                        0.0
        BC_eta_1                      7.8539816339744828e-01
        BC_mean_x_1                   0.0
        BC_sigma_x_1                  {source_sigma:.16e}
        BC_G_1                        {source_rate:.16e}
        BC_G_fileName_1               {tag}_one.txt
        BC_G_NS_1                     80

        SPECIES2                      1
        NPC2                          {npc}
        pctSupPartOutput2             100.0
        Z2                            -1
        M2                            5.4857990904410000e-04

        IC_type_2                     1
        IC_Tper_2                     5.0000000000000000e+03
        IC_Tper_fileName_2            {tag}_one.txt
        IC_Tper_NX_2                  80
        IC_Tpar_2                     5.0000000000000000e+03
        IC_Tpar_fileName_2            {tag}_one.txt
        IC_Tpar_NX_2                  80
        IC_densityFraction_2          1.0
        IC_densityFraction_fileName_2 {tag}_one.txt
        IC_densityFraction_NX_2       80
        IC_weightScale_2              0.0

        BC_type_2                     1
        BC_T_2                        5.0000000000000000e+03
        BC_E_2                        0.0
        BC_eta_2                      7.8539816339744828e-01
        BC_mean_x_2                   0.0
        BC_sigma_x_2                  {source_sigma:.16e}
        BC_G_2                        {source_rate:.16e}
        BC_G_fileName_2               {tag}_one.txt
        BC_G_NS_2                     80
        """
    ).strip() + "\n"


def numeric_keys(handle: h5py.File) -> list[str]:
    return sorted([key for key in handle.keys() if key.isdigit()], key=lambda item: int(item))


def read_scalar(handle: h5py.File, dataset: str) -> float:
    return float(np.asarray(handle[dataset]).reshape(-1)[0])


def read_vector(handle: h5py.File, dataset: str) -> np.ndarray | None:
    if dataset not in handle:
        return None
    return np.asarray(handle[dataset], dtype=float).reshape(-1)


def prepare_case(root: Path, case: SourceCase, npc: int, source_rate: float, source_sigma: float, sim_gyroperiods: float, output_gyroperiods: float) -> Path:
    picos_files = root / case.tag / "picosFILES"
    input_dir = picos_files / "inputFiles"
    input_dir.mkdir(parents=True, exist_ok=True)
    deck_source_rate = source_rate if case.pair_source else 0.0
    (input_dir / f"input_file_{case.tag}.input").write_text(input_deck(case.tag, case.pair_source, deck_source_rate, sim_gyroperiods, output_gyroperiods))
    (input_dir / f"ions_properties_{case.tag}.ion").write_text(ion_deck(case.tag, npc, source_rate, source_sigma))
    write_vector(input_dir / f"{case.tag}_one.txt", np.ones(80))
    return picos_files


def run_case(binary: Path, picos_files: Path, case: SourceCase, mpi_ranks: int, timeout: int) -> None:
    output_root = picos_files / "outputFiles"
    if output_root.exists():
        shutil.rmtree(output_root)
    command = ["mpirun", "-np", str(mpi_ranks), str(binary), "1-D", "outputFiles", case.tag]
    env = os.environ.copy()
    env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    log_path = picos_files.parent / "run.log"
    with log_path.open("w") as stream:
        subprocess.run(command, cwd=picos_files, env=env, stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=timeout)


def read_case(hdf5_dir: Path) -> dict[str, object]:
    main_path = hdf5_dir / "main.h5"
    particles_path = hdf5_dir / "PARTICLES_FILE_0.h5"
    if not main_path.exists():
        raise FileNotFoundError(main_path)
    if not particles_path.exists():
        raise FileNotFoundError(particles_path)

    with h5py.File(main_path, "r") as handle:
        z = np.asarray(handle["geometry/x_m"], dtype=float).reshape(-1)

    data: dict[str, object] = {"z": z, "species": {}}
    with h5py.File(particles_path, "r") as handle:
        steps = numeric_keys(handle)
        times = np.asarray([read_scalar(handle, f"{step}/time") for step in steps], dtype=float)
        data["times"] = times
        species_group = handle[f"{steps[0]}/ions"]
        species_names = sorted(species_group.keys(), key=lambda item: int(item.split("_")[-1]))
        for species in species_names:
            species_data: dict[str, np.ndarray] = {}
            for variable in ("n_m", "Tpar_m", "Tper_m"):
                columns = []
                for step in steps:
                    values = read_vector(handle, f"{step}/ions/{species}/{variable}")
                    if values is None:
                        values = read_vector(handle, f"{step}/ions/{species}/{variable}/x")
                    if values is None:
                        columns = []
                        break
                    columns.append(values)
                if columns:
                    species_data[variable] = np.column_stack(columns)
            data["species"][species] = species_data

        boundary_n5 = []
        for step in steps:
            n5 = read_vector(handle, f"{step}/boundary/N5")
            boundary_n5.append(float(n5[0]) if n5 is not None and n5.size else np.nan)
        data["N5"] = np.asarray(boundary_n5, dtype=float)
    return data


def trim_to_common_grid(data: dict[str, object]) -> None:
    z = data["z"]
    species_data = data["species"]
    n = min([len(z)] + [arr.shape[0] for per_species in species_data.values() for arr in per_species.values()])
    data["z"] = z[:n]
    for per_species in species_data.values():
        for key, values in list(per_species.items()):
            per_species[key] = values[:n, :]


def total_particles(z: np.ndarray, density: np.ndarray, radius: float = 0.02) -> np.ndarray:
    area = np.pi * radius * radius
    return np.asarray([np.trapz(density[:, ii], z) * area for ii in range(density.shape[1])])


def coverage_fraction(density: np.ndarray, threshold_fraction: float = 1.0e-3) -> float:
    final = density[:, -1]
    peak = float(np.nanmax(final)) if final.size else 0.0
    if peak <= 0.0:
        return 0.0
    return float(np.count_nonzero(final > threshold_fraction * peak) / final.size)


def finite_temperature_median(density: np.ndarray, temperature: np.ndarray) -> float:
    final_density = density[:, -1]
    final_temperature = temperature[:, -1]
    peak = float(np.nanmax(final_density)) if final_density.size else 0.0
    mask = np.isfinite(final_temperature) & (final_density > max(peak * 1.0e-3, 0.0))
    if not np.any(mask):
        return float("nan")
    return float(np.nanmedian(final_temperature[mask]))


def summarize_case(case: SourceCase, data: dict[str, object]) -> list[dict[str, object]]:
    z = data["z"]
    times = data["times"]
    rows: list[dict[str, object]] = []
    for species, species_data in data["species"].items():
        density = species_data["n_m"]
        n_total = total_particles(z, density)
        tpar = species_data.get("Tpar_m")
        tper = species_data.get("Tper_m")
        peak_index = int(np.nanargmax(density[:, -1])) if density.size else 0
        rows.append(
            {
                "case": case.tag,
                "source_model": "pair_explicit_rate" if case.pair_source else "legacy",
                "species": SPECIES_LABEL.get(species, species),
                "t_initial_s": float(times[0]),
                "t_final_s": float(times[-1]),
                "N_initial": float(n_total[0]),
                "N_final": float(n_total[-1]),
                "N_growth": float(n_total[-1] - n_total[0]),
                "density_peak_final_m-3": float(np.nanmax(density[:, -1])),
                "density_peak_z_m": float(z[peak_index]),
                "density_coverage_fraction": coverage_fraction(density),
                "Tpar_median_final_eV": finite_temperature_median(density, tpar) if tpar is not None else float("nan"),
                "Tper_median_final_eV": finite_temperature_median(density, tper) if tper is not None else float("nan"),
            }
        )
    return rows


def plot_case(case: SourceCase, data: dict[str, object], out_dir: Path) -> None:
    z = data["z"]
    times_us = data["times"] * 1.0e6
    species_data = data["species"]
    species_names = sorted(species_data.keys(), key=lambda item: int(item.split("_")[-1]))

    fig, axes = plt.subplots(len(species_names), 3, figsize=(18, 5.0 * len(species_names)), sharex=True)
    if len(species_names) == 1:
        axes = np.asarray([axes])
    for row, species in enumerate(species_names):
        label = SPECIES_LABEL.get(species, species)
        for col, variable in enumerate(("n_m", "Tpar_m", "Tper_m")):
            ax = axes[row, col]
            values = species_data[species][variable]
            plot_values = values.T
            if variable == "n_m":
                vmax = float(np.nanpercentile(plot_values, 99.5)) if np.any(np.isfinite(plot_values)) else 1.0
                im = ax.pcolormesh(z, times_us, plot_values, shading="auto", cmap="viridis", vmin=0.0, vmax=vmax)
                cbar_label = "n [m^-3]"
            else:
                masked = np.where(species_data[species]["n_m"].T > np.nanmax(species_data[species]["n_m"]) * 1.0e-4, plot_values, np.nan)
                im = ax.pcolormesh(z, times_us, masked, shading="auto", cmap="plasma", vmin=0.0, vmax=7000.0)
                cbar_label = "T [eV]"
            ax.axvline(0.0, color="w" if variable == "n_m" else "k", linestyle="--", linewidth=1.2)
            ax.set_title(f"{label} {variable}")
            ax.set_ylabel("time [us]")
            ax.grid(True, alpha=0.18)
            cbar = fig.colorbar(im, ax=ax, pad=0.01)
            cbar.set_label(cbar_label)
    for ax in axes[-1, :]:
        ax.set_xlabel("z [m]")
    fig.suptitle(f"Source build-up validation: {case.tag}", fontsize=18)
    fig.tight_layout()
    fig.savefig(out_dir / f"{case.tag}_density_temperature_build_up.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex="col")
    sample_indices = sorted(set([0, max(0, len(times_us) // 4), max(0, len(times_us) // 2), len(times_us) - 1]))
    for species in species_names:
        label = SPECIES_LABEL.get(species, species)
        density = species_data[species]["n_m"]
        n_total = total_particles(z, density)
        axes[0, 0].plot(times_us, n_total, linewidth=2.0, label=label)
        for idx in sample_indices:
            axes[0, 1].plot(z, density[:, idx], linewidth=1.8, label=f"{label}, {times_us[idx]:.2f} us")
        axes[1, 0].plot(times_us, np.nanmax(density, axis=0), linewidth=2.0, label=label)
        axes[1, 1].plot(times_us, [coverage_fraction(density[:, :ii + 1]) for ii in range(density.shape[1])], linewidth=2.0, label=label)
    axes[0, 0].set_ylabel("N in domain")
    axes[0, 1].set_ylabel("n [m^-3]")
    axes[1, 0].set_ylabel("peak n [m^-3]")
    axes[1, 1].set_ylabel("coverage fraction")
    axes[1, 0].set_xlabel("time [us]")
    axes[1, 1].set_xlabel("time [us]")
    axes[0, 1].set_xlabel("z [m]")
    axes[0, 0].grid(True, alpha=0.25)
    axes[0, 1].grid(True, alpha=0.25)
    axes[1, 0].grid(True, alpha=0.25)
    axes[1, 1].grid(True, alpha=0.25)
    for ax in axes.reshape(-1):
        ax.legend(fontsize=8)
    fig.suptitle(f"Source build-up metrics: {case.tag}", fontsize=18)
    fig.tight_layout()
    fig.savefig(out_dir / f"{case.tag}_source_build_up_metrics.png", dpi=220)
    plt.close(fig)


def write_summary(output_dir: Path, rows: list[dict[str, object]], figures: list[Path]) -> None:
    csv_path = output_dir / "source_build_up_summary.csv"
    fieldnames = sorted({key for row in rows for key in row})
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# PICOS++ source build-up validation",
        "",
        "The source validation starts computational markers with zero physical weight (`IC_weightScale=0`) and uses the z=0 warm source to populate the domain.",
        "",
        "## Metrics",
        "",
    ]
    for row in rows:
        lines.append(
            f"- {row['case']} {row['species']}: N_final={row['N_final']:.6e}, "
            f"peak n={row['density_peak_final_m-3']:.6e} m^-3 at z={row['density_peak_z_m']:.4f} m, "
            f"coverage={row['density_coverage_fraction']:.3f}, "
            f"Tpar={row['Tpar_median_final_eV']:.2f} eV, Tper={row['Tper_median_final_eV']:.2f} eV"
        )
    lines.extend(["", "## Figures", ""])
    for figure in figures:
        lines.append(f"- `{figure}`")
    lines.append("")
    (output_dir / "source_build_up_summary.md").write_text("\n".join(lines))


def validate_summary(rows: list[dict[str, object]], min_coverage: float) -> None:
    failures: list[str] = []
    for row in rows:
        if row["N_growth"] <= 0.0:
            failures.append(f"{row['case']} {row['species']} did not increase total particle content")
        if row["density_peak_final_m-3"] <= 0.0:
            failures.append(f"{row['case']} {row['species']} has zero final peak density")
        if row["density_coverage_fraction"] < min_coverage:
            failures.append(
                f"{row['case']} {row['species']} coverage {row['density_coverage_fraction']:.3f} "
                f"is below required {min_coverage:.3f}"
            )
        if not np.isfinite(row["Tpar_median_final_eV"]) or not np.isfinite(row["Tper_median_final_eV"]):
            failures.append(f"{row['case']} {row['species']} has non-finite final temperature moment")
    if failures:
        raise RuntimeError("Source build-up validation failed:\n  - " + "\n  - ".join(failures))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=PICOS_ROOT / "build" / "picosFILES" / "src" / "xpicos")
    parser.add_argument("--mpi-ranks", type=int, default=2)
    parser.add_argument("--npc", type=int, default=8)
    parser.add_argument("--source-rate", type=float, default=1.0e22)
    parser.add_argument("--source-sigma", type=float, default=0.04)
    parser.add_argument("--sim-gyroperiods", type=float, default=0.5)
    parser.add_argument("--output-gyroperiods", type=float, default=0.05)
    parser.add_argument("--min-coverage", type=float, default=0.8)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--skip-runs", action="store_true")
    parser.add_argument("--case", choices=["legacy", "pair", "both"], default="both")
    args = parser.parse_args()

    binary = args.binary.resolve()
    if not binary.exists():
        raise FileNotFoundError(binary)
    if args.mpi_ranks % 2 != 0:
        raise ValueError("PICOS++ requires an even number of MPI ranks.")

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    output_dir = PICOS_ROOT / "validation" / f"source_build_up_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)

    cases: list[SourceCase] = []
    if args.case in ("legacy", "both"):
        cases.append(SourceCase("source_build_up_legacy_z0_empty", False))
    if args.case in ("pair", "both"):
        cases.append(SourceCase("source_build_up_pair_explicit_rate_z0_empty", True))

    rows: list[dict[str, object]] = []
    figures: list[Path] = []
    for case in cases:
        picos_files = prepare_case(output_dir, case, args.npc, args.source_rate, args.source_sigma, args.sim_gyroperiods, args.output_gyroperiods)
        if not args.skip_runs:
            print(f"running {case.tag}")
            run_case(binary, picos_files, case, args.mpi_ranks, args.timeout)
        hdf5_dir = picos_files / "outputFiles" / case.tag / "HDF5"
        data = read_case(hdf5_dir)
        trim_to_common_grid(data)
        rows.extend(summarize_case(case, data))
        plot_case(case, data, output_dir)
        figures.append(output_dir / f"{case.tag}_density_temperature_build_up.png")
        figures.append(output_dir / f"{case.tag}_source_build_up_metrics.png")

    write_summary(output_dir, rows, figures)
    validate_summary(rows, args.min_coverage)
    print(output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
