#!/usr/bin/env python3
"""Compare the PICOS_ECH particle source against origin/master.

The test is intentionally small and artificial: particles stream in a short
1D domain with no collisions, RF, or field solve.  Boundary-loss reinjection is
forced often enough that the source location can be inferred from jumps in
saved particle positions.  The default source is a narrow Gaussian at z=0.
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
import matplotlib.pyplot as plt
import numpy as np


PICOS_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MASTER_ROOT = PICOS_ROOT.parent / "PICOS_master_source_compare"


@dataclass(frozen=True)
class Case:
    name: str
    branch_label: str
    exe: Path
    pair_source: bool
    pair_mode: int
    pair_profile_file: str


def write_vector(path: Path, values: np.ndarray) -> None:
    np.savetxt(path, np.asarray(values, dtype=float), fmt="%.16e")


def input_deck(case: Case) -> str:
    pair_switch = 1 if case.pair_source else 0
    return textwrap.dedent(
        f"""
        // PICOS source-location validation deck.
        // Purpose: compare source injection against origin/master.
        // Physics: no RF, no collisions, no E-field solve; particles stream and
        // reinject from a narrow source centered at z=0.
        // =============================================================================
        mpisForFields               1
        quietStart                  1
        IC_velocityDistributionModel 0
        IC_randomSeed               137
        numberOfParticleSpecies     2
        numberOfTracerSpecies       0
        advanceParticleMethod       1

        // Characteristic values:
        // =============================================================================
        CV_ne                       1.0000000000000000e+17
        CV_Te                       5.0000000000000000e+03
        CV_B                        1.0000000000000001e-01
        CV_Tpar                     5.0000000000000000e+03
        CV_Tper                     5.0000000000000000e+03

        // Simulation time is in background-ion gyroperiod units.
        // =============================================================================
        DTc                         1.0000000000000001e-01
        simulationTime              2.0000000000000000e+00

        // Switches:
        // =============================================================================
        SW_EfieldSolve              0
        SW_fieldSolveModel          0
        SW_BfieldSolve              0
        SW_Collisions               0
        CollOperType                2
        SW_collisionConservationProjection 0
        collisionRandomSeed         2468
        SW_RFheating                0
        SW_RFheatingIons            0
        SW_RFheatingElectrons       0
        SW_pairSource               {pair_switch}
        SW_relativisticElectrons    0
        SW_electronGyroTimeStepLimiter 0
        SW_electronPlasmaTimeStepLimiter 0
        SW_advancePos               1
        SW_linearSolve              0

        // Magnetic field initial conditions:
        // =============================================================================
        IC_uniformBfield            1
        IC_BX                       1.0000000000000001e-01
        IC_BY                       0.0
        IC_BZ                       0.0
        IC_BX_NX                    34
        IC_BX_fileName              source_compare_unit_profile.txt
        IC_phiLeft                  0.0
        IC_phiRight                 0.0
        Poisson_BCModel             1
        Poisson_sheathCoefficient   3.0
        ReformulatedPoisson_lambda  -1.0
        ReformulatedPoisson_quasiNeutral 0

        // Geometry:
        // =============================================================================
        dp                          1.0000000000000001e-01
        r1                          0.0
        r2                          2.0000000000000000e-02
        LX_min                      -5.0000000000000000e-01
        LX_max                      5.0000000000000000e-01

        // Electron fluid/profile initial conditions:
        // =============================================================================
        IC_ne                       1.0000000000000000e+17
        IC_Te                       5.0000000000000000e+03
        IC_Te_NX                    34
        IC_Te_fileName              source_compare_unit_profile.txt

        // Coupled electron-ion pair source.  origin/master ignores these keys.
        // positionMode: 0 = Gaussian, 1 = profile-weighted.
        // =============================================================================
        pairSource_ionSpecies       1
        pairSource_electronSpecies  2
        pairSource_rate             1.0000000000000000e+22
        pairSource_mean_x           0.0
        pairSource_sigma_x          2.0000000000000000e-02
        pairSource_Ti_birth         5.0000000000000000e+03
        pairSource_Te_birth         5.0000000000000000e+03
        pairSource_Ei_birth         0.0
        pairSource_Ee_birth         0.0
        pairSource_eta_i            0.7853981633974483
        pairSource_eta_e            0.7853981633974483
        pairSource_positionMode     {case.pair_mode}
        pairSource_fileName         {case.pair_profile_file}
        pairSource_NS               34
        pairSource_maxParticleWeight 1000

        // Legacy RF keys required by origin/master, disabled here.
        // =============================================================================
        RF_Prf                      0.0
        RF_n_harmonic               1
        RF_freq                     1.0
        RF_x1                       -5.0000000000000000e-01
        RF_x2                       5.0000000000000000e-01
        RF_t_ON                     0.0
        RF_t_OFF                    2.0000000000000000e+00
        RF_kpar                     0.0
        RF_kper                     0.0
        RF_handedness               -1
        RF_Prf_fileName             source_compare_unit_profile.txt
        RF_Prf_NS                   34

        // Species-split RF keys used by PICOS_ECH, disabled here.
        // =============================================================================
        RF_ion_Prf                      0.0
        RF_ion_n_harmonic               1
        RF_ion_freq                     1.0
        RF_ion_x1                       -5.0000000000000000e-01
        RF_ion_x2                       5.0000000000000000e-01
        RF_ion_t_ON                     0.0
        RF_ion_t_OFF                    2.0000000000000000e+00
        RF_ion_kpar                     0.0
        RF_ion_kper                     0.0
        RF_ion_handedness               -1
        RF_ion_EfieldMode               1
        RF_ion_EfieldAmplitude          0.0
        RF_ion_maxEnergyGainFraction    0.0
        RF_ion_maxParticleEnergy        0.0
        RF_ion_maxVelocityFractionC     0.0
        RF_ion_Prf_fileName             source_compare_unit_profile.txt
        RF_ion_Prf_NS                   34

        RF_electron_Prf                      0.0
        RF_electron_n_harmonic               1
        RF_electron_freq                     1.0
        RF_electron_x1                       -5.0000000000000000e-01
        RF_electron_x2                       5.0000000000000000e-01
        RF_electron_t_ON                     0.0
        RF_electron_t_OFF                    2.0000000000000000e+00
        RF_electron_kpar                     0.0
        RF_electron_kper                     0.0
        RF_electron_handedness               -1
        RF_electron_EfieldMode               1
        RF_electron_EfieldAmplitude          0.0
        RF_electron_maxEnergyGainFraction    0.0
        RF_electron_maxParticleEnergy        0.0
        RF_electron_maxVelocityFractionC     0.0
        RF_electron_Prf_fileName             source_compare_unit_profile.txt
        RF_electron_Prf_NS                   34

        // Output variables:
        // =============================================================================
        outputCadence               1.0000000000000000e-03
        outputs_variables           {{X_p,V_p,a_p,BX_p,n_m,Tpar_m,Tper_m,u_m}}

        // Data smoothing:
        // =============================================================================
        smoothingParameter          0.0
        filtersPerIterationFields   0
        filtersPerIterationIons     0
        """
    ).strip() + "\n"


def ion_deck() -> str:
    return textwrap.dedent(
        """
        // PICOS source-location validation species deck.
        // Species 1 is D+.  Species 2 is kinetic electron-like.
        // =============================================================================
        SPECIES1                      1
        NPC1                          32
        pctSupPartOutput1             100.0
        Z1                            1
        M1                            2.0000000000000000e+00

        IC_type_1                     1
        IC_Tper_1                     5.0000000000000000e+03
        IC_Tper_fileName_1            source_compare_unit_profile.txt
        IC_Tper_NX_1                  34
        IC_Tpar_1                     5.0000000000000000e+03
        IC_Tpar_fileName_1            source_compare_unit_profile.txt
        IC_Tpar_NX_1                  34
        IC_densityFraction_1          1.0
        IC_densityFraction_fileName_1 source_compare_unit_profile.txt
        IC_densityFraction_NX_1       34

        BC_type_1                     1
        BC_T_1                        5.0000000000000000e+03
        BC_E_1                        0.0
        BC_eta_1                      0.7853981633974483
        BC_mean_x_1                   0.0
        BC_sigma_x_1                  2.0000000000000000e-02
        BC_G_1                        1.0000000000000000e+22
        BC_G_fileName_1               source_compare_unit_profile.txt
        BC_G_NS_1                     34

        SPECIES2                      1
        NPC2                          32
        pctSupPartOutput2             100.0
        Z2                            -1
        M2                            5.4857990904410000e-04

        IC_type_2                     1
        IC_Tper_2                     5.0000000000000000e+03
        IC_Tper_fileName_2            source_compare_unit_profile.txt
        IC_Tper_NX_2                  34
        IC_Tpar_2                     5.0000000000000000e+03
        IC_Tpar_fileName_2            source_compare_unit_profile.txt
        IC_Tpar_NX_2                  34
        IC_densityFraction_2          1.0
        IC_densityFraction_fileName_2 source_compare_unit_profile.txt
        IC_densityFraction_NX_2       34

        BC_type_2                     1
        BC_T_2                        5.0000000000000000e+03
        BC_E_2                        0.0
        BC_eta_2                      0.7853981633974483
        BC_mean_x_2                   0.0
        BC_sigma_x_2                  2.0000000000000000e-02
        BC_G_2                        1.0000000000000000e+22
        BC_G_fileName_2               source_compare_unit_profile.txt
        BC_G_NS_2                     34
        """
    ).strip() + "\n"


def prepare_case(case: Case, base_dir: Path) -> Path:
    picos_files = base_dir / case.name / "picosFILES"
    input_dir = picos_files / "inputFiles"
    input_dir.mkdir(parents=True, exist_ok=True)

    (input_dir / "input_file.input").write_text(input_deck(case))
    (input_dir / "ions_properties.ion").write_text(ion_deck())
    write_vector(input_dir / "source_compare_unit_profile.txt", np.ones(34))

    z = np.linspace(-0.5, 0.5, 34)
    source_profile = np.exp(-0.5 * (z / 0.02) ** 2)
    write_vector(input_dir / "source_compare_pair_profile_z0.txt", source_profile)

    return picos_files


def run_case(case: Case, picos_files: Path, timeout: int) -> None:
    out_dir = picos_files / "outputFiles"
    cmd = [
        "mpirun",
        "-np",
        "2",
        str(case.exe),
        "source_compare",
        str(out_dir),
    ]
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    log = picos_files.parent / "run.log"
    with log.open("w") as stream:
        subprocess.run(
            cmd,
            cwd=picos_files,
            env=env,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=True,
            timeout=timeout,
        )


def sorted_snapshot_names(h5file: h5py.File) -> list[str]:
    names = [key for key in h5file.keys() if key.isdigit()]
    return sorted(names, key=lambda item: int(item))


def find_dataset_paths(h5file: h5py.File, suffix: str) -> list[str]:
    matches: list[str] = []

    def visitor(name: str, obj: h5py.Dataset) -> None:
        if isinstance(obj, h5py.Dataset) and name.endswith(suffix):
            matches.append(name)

    h5file.visititems(visitor)
    return matches


def read_birth_jumps(particle_file: Path, species_name: str) -> np.ndarray:
    with h5py.File(particle_file, "r") as handle:
        snapshots = sorted_snapshot_names(handle)
        previous = None
        births: list[np.ndarray] = []
        for snap in snapshots:
            path = f"{snap}/ions/{species_name}/X_p"
            if path not in handle:
                continue
            x = np.asarray(handle[path][...], dtype=float).reshape(-1)
            if previous is not None and previous.shape == x.shape:
                dx = x - previous
                jumps = np.abs(dx) > 0.15
                if np.any(jumps):
                    births.append(x[jumps])
            previous = x
    if births:
        return np.concatenate(births)
    return np.empty(0)


def read_n5_stats(hdf_dir: Path) -> dict[str, float | int]:
    values: list[float] = []
    for particle_file in sorted(hdf_dir.glob("PARTICLES_FILE_*.h5")):
        with h5py.File(particle_file, "r") as handle:
            paths = find_dataset_paths(handle, "boundary/N5")
            for path in paths:
                value = float(np.asarray(handle[path][...], dtype=float).reshape(-1)[0])
                if np.isfinite(value):
                    values.append(value)
    if not values:
        return {
            "n5_output_count": 0,
            "n5_nonzero_count": 0,
            "n5_sum_rate_samples": float("nan"),
            "n5_mean_nonzero_rate_s-1": float("nan"),
            "n5_max_rate_s-1": float("nan"),
        }

    arr = np.asarray(values)
    nonzero = arr[np.abs(arr) > 0.0]
    return {
        "n5_output_count": int(arr.size),
        "n5_nonzero_count": int(nonzero.size),
        "n5_sum_rate_samples": float(np.sum(arr)),
        "n5_mean_nonzero_rate_s-1": float(np.mean(nonzero)) if nonzero.size else 0.0,
        "n5_max_rate_s-1": float(np.max(arr)),
    }


def read_final_density_peak(hdf_dir: Path, species_name: str) -> tuple[float, float]:
    particle_file = hdf_dir / "PARTICLES_FILE_0.h5"
    field_file = hdf_dir / "FIELDS_FILE_0.h5"
    if not particle_file.exists():
        return float("nan"), float("nan")

    with h5py.File(particle_file, "r") as particles:
        snapshots = sorted_snapshot_names(particles)
        if not snapshots:
            return float("nan"), float("nan")
        snap = snapshots[-1]
        density_path = f"{snap}/ions/{species_name}/n_m/x"
        if density_path not in particles:
            density_path = f"{snap}/ions/{species_name}/n_m"
        if density_path not in particles:
            return float("nan"), float("nan")
        density = np.asarray(particles[density_path][...], dtype=float).reshape(-1)

    x = np.linspace(-0.5, 0.5, density.size)
    if field_file.exists():
        with h5py.File(field_file, "r") as fields:
            snapshots = sorted_snapshot_names(fields)
            if snapshots:
                x_path = f"{snapshots[-1]}/fields/BX_m/x"
                if x_path in fields:
                    bx = np.asarray(fields[x_path][...], dtype=float).reshape(-1)
                    if bx.size == density.size:
                        x = np.linspace(-0.5, 0.5, bx.size)

    peak_index = int(np.nanargmax(density))
    return float(x[peak_index]), float(density[peak_index])


def summarize_case(case: Case, hdf_dir: Path) -> dict[str, float | str]:
    particle_file = hdf_dir / "PARTICLES_FILE_0.h5"
    ion_births = read_birth_jumps(particle_file, "spp_1") if particle_file.exists() else np.empty(0)
    electron_births = read_birth_jumps(particle_file, "spp_2") if particle_file.exists() else np.empty(0)
    ion_peak_z, ion_peak_n = read_final_density_peak(hdf_dir, "spp_1")
    electron_peak_z, electron_peak_n = read_final_density_peak(hdf_dir, "spp_2")

    def stats(prefix: str, values: np.ndarray) -> dict[str, float | int]:
        if values.size == 0:
            return {
                f"{prefix}_birth_count": 0,
                f"{prefix}_birth_mean_z_m": float("nan"),
                f"{prefix}_birth_std_z_m": float("nan"),
                f"{prefix}_birth_abs_p95_m": float("nan"),
            }
        return {
            f"{prefix}_birth_count": int(values.size),
            f"{prefix}_birth_mean_z_m": float(np.mean(values)),
            f"{prefix}_birth_std_z_m": float(np.std(values)),
            f"{prefix}_birth_abs_p95_m": float(np.percentile(np.abs(values), 95)),
        }

    row: dict[str, float | str] = {
        "case": case.name,
        "branch": case.branch_label,
        "source_model": "pair" if case.pair_source else "legacy",
        "pair_position_mode": case.pair_mode if case.pair_source else -1,
        "ion_density_peak_z_m": ion_peak_z,
        "ion_density_peak_m3": ion_peak_n,
        "electron_density_peak_z_m": electron_peak_z,
        "electron_density_peak_m3": electron_peak_n,
    }
    row.update(read_n5_stats(hdf_dir))
    row.update(stats("ion", ion_births))
    row.update(stats("electron", electron_births))
    return row


def plot_birth_distributions(rows: list[dict[str, float | str]], run_dirs: dict[str, Path], output_dir: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
    bins = np.linspace(-0.12, 0.12, 49)

    for row in rows:
        name = str(row["case"])
        hdf_dir = run_dirs[name] / "outputFiles" / "HDF5"
        particle_file = hdf_dir / "PARTICLES_FILE_0.h5"
        if not particle_file.exists():
            continue
        for ax, species, label_prefix in [
            (axes[0], "spp_1", "D+"),
            (axes[1], "spp_2", "e-"),
        ]:
            births = read_birth_jumps(particle_file, species)
            if births.size:
                ax.hist(
                    births,
                    bins=bins,
                    histtype="step",
                    density=True,
                    linewidth=2.0,
                    label=name,
                )
            ax.set_ylabel(f"{label_prefix} birth PDF")
            ax.axvline(0.0, color="k", linestyle="--", linewidth=1.0)
            ax.grid(True, alpha=0.25)
    axes[-1].set_xlabel("inferred birth z [m]")
    for ax in axes:
        ax.legend(fontsize=8)
    fig.suptitle("PICOS source-position comparison, z0 source decks")
    fig.tight_layout()
    fig.savefig(output_dir / "source_birth_position_compare.png", dpi=220)
    plt.close(fig)


def analyze_mpex_source_profile(picos_root: Path, output_dir: Path) -> dict[str, float | str]:
    profile_name = "mpex_scenario14_ex8_steady_2ms_p262144_coll_mpexprof_nersc_nonrel_pair_source_norm.txt"
    profile_path = picos_root / "picosFILES" / "inputFiles" / profile_name
    if not profile_path.exists():
        return {"mpex_profile_file": str(profile_path), "mpex_profile_exists": 0}
    values = np.loadtxt(profile_path, dtype=float)
    z = np.linspace(-2.0, 8.0, values.size)
    positive = np.clip(values, 0.0, None)
    weighted_mean = float(np.sum(z * positive) / np.sum(positive))
    peak_z = float(z[int(np.argmax(positive))])

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(z, positive, color="tab:red", linewidth=2.2)
    ax.axvline(0.0, color="k", linestyle="--", linewidth=1.2, label="z=0")
    ax.axvline(peak_z, color="tab:red", linestyle=":", linewidth=1.8, label=f"peak z={peak_z:.3f} m")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("source profile [arb.]")
    ax.set_title("Current MPEX pair-source profile in PICOS_ECH input")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "current_mpex_pair_source_profile.png", dpi=220)
    plt.close(fig)

    return {
        "mpex_profile_file": str(profile_path),
        "mpex_profile_exists": 1,
        "mpex_profile_peak_z_m": peak_z,
        "mpex_profile_weighted_mean_z_m": weighted_mean,
        "mpex_profile_n_points": int(values.size),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current-root", type=Path, default=PICOS_ROOT)
    parser.add_argument("--master-root", type=Path, default=DEFAULT_MASTER_ROOT)
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--skip-runs", action="store_true")
    args = parser.parse_args()

    current_exe = args.current_root / "build_nopull" / "picosFILES" / "src" / "xpicos"
    master_exe = args.master_root / "build_nopull" / "picosFILES" / "src" / "xpicos"
    for exe in [current_exe, master_exe]:
        if not exe.exists():
            raise FileNotFoundError(f"Missing executable: {exe}")

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    output_dir = args.current_root / "validation" / f"source_compare_z0_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        Case("master_legacy_z0", "origin/master", master_exe, False, 0, "source_compare_unit_profile.txt"),
        Case("current_legacy_z0", "PICOS_ECH", current_exe, False, 0, "source_compare_unit_profile.txt"),
        Case("current_pair_gaussian_z0", "PICOS_ECH", current_exe, True, 0, "source_compare_unit_profile.txt"),
        Case("current_pair_profile_z0", "PICOS_ECH", current_exe, True, 1, "source_compare_pair_profile_z0.txt"),
    ]

    run_dirs: dict[str, Path] = {}
    for case in cases:
        picos_files = prepare_case(case, output_dir)
        run_dirs[case.name] = picos_files
        if not args.skip_runs:
            print(f"running {case.name}")
            run_case(case, picos_files, args.timeout)

    rows = []
    for case in cases:
        hdf_dir = run_dirs[case.name] / "outputFiles" / "HDF5"
        rows.append(summarize_case(case, hdf_dir))

    mpex_profile_summary = analyze_mpex_source_profile(args.current_root, output_dir)

    csv_path = output_dir / "source_compare_summary.csv"
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with (output_dir / "mpex_pair_source_profile_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(mpex_profile_summary.keys()))
        writer.writeheader()
        writer.writerow(mpex_profile_summary)

    plot_birth_distributions(rows, run_dirs, output_dir)

    print(f"wrote {csv_path}")
    print(f"wrote {output_dir / 'source_birth_position_compare.png'}")
    print(f"wrote {output_dir / 'current_mpex_pair_source_profile.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
