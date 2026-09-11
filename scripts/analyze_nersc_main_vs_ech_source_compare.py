#!/usr/bin/env python3
"""Compare NERSC hybrid source runs from PICOS main and PICOS_ECH branches."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_RUN_ROOT = Path("PICOS_NERSC_main_vs_ech_hybrid_source_1ms")
CASES = {
    "main": {
        "subdir": "main_hybrid_1ms",
        "tag": "picos_main_hybrid_source_1ms",
        "label": "PICOS main",
    },
    "ech": {
        "subdir": "ech_branch_hybrid_1ms",
        "tag": "picos_ech_branch_hybrid_source_1ms",
        "label": "PICOS_ECH",
    },
}


def numeric_keys(handle: h5py.File) -> list[str]:
    return sorted([key for key in handle.keys() if key.isdigit()], key=lambda item: int(item))


def first_dataset(handle: h5py.File, candidates: list[str]) -> np.ndarray | None:
    for path in candidates:
        if path in handle:
            return np.asarray(handle[path], dtype=float).reshape(-1)
    return None


def hdf_dir(run_root: Path, case: dict[str, str]) -> Path:
    path = (
        run_root
        / case["subdir"]
        / "run"
        / "picosFILES"
        / "outputFiles"
        / case["tag"]
        / "HDF5"
    )
    if not path.is_dir():
        raise FileNotFoundError(f"Missing HDF5 output directory: {path}")
    return path


def read_case(path: Path) -> dict[str, np.ndarray]:
    main_path = path / "main.h5"
    particle_path = path / "PARTICLES_FILE_0.h5"
    field_path = path / "FIELDS_FILE_0.h5"
    if not main_path.exists():
        raise FileNotFoundError(main_path)
    if not particle_path.exists():
        raise FileNotFoundError(particle_path)

    with h5py.File(main_path, "r") as handle:
        z = first_dataset(handle, ["geometry/x_m", "x_m"])
        if z is None:
            raise KeyError(f"No geometry/x_m dataset in {main_path}")

    out: dict[str, np.ndarray] = {"z": z}
    with h5py.File(particle_path, "r") as handle:
        snaps = numeric_keys(handle)
        if not snaps:
            raise KeyError(f"No numeric snapshots in {particle_path}")
        times = np.asarray([float(np.asarray(handle[f"{snap}/time"]).reshape(-1)[0]) for snap in snaps])
        out["time"] = times
        for variable in ("n_m", "Tpar_m", "Tper_m", "u_m"):
            columns = []
            for snap in snaps:
                values = first_dataset(
                    handle,
                    [
                        f"{snap}/ions/spp_1/{variable}",
                        f"{snap}/ions/spp_1/{variable}/x",
                    ],
                )
                if values is None:
                    columns = []
                    break
                columns.append(values)
            if columns:
                out[variable] = np.column_stack(columns)

    if field_path.exists():
        with h5py.File(field_path, "r") as handle:
            snaps = numeric_keys(handle)
            for variable in ("BX_m", "EX_m"):
                columns = []
                for snap in snaps:
                    values = first_dataset(
                        handle,
                        [
                            f"{snap}/fields/{variable}",
                            f"{snap}/fields/{variable}/x",
                        ],
                    )
                    if values is None:
                        columns = []
                        break
                    columns.append(values)
                if columns:
                    out[variable] = np.column_stack(columns)
    return out


def trim_common(data: dict[str, dict[str, np.ndarray]]) -> None:
    nz = min(len(case["z"]) for case in data.values())
    nt = min(case["time"].size for case in data.values())
    for case in data.values():
        case["z"] = case["z"][:nz]
        case["time"] = case["time"][:nt]
        for key, value in list(case.items()):
            if value.ndim == 2:
                case[key] = value[:nz, :nt]


def rel_l2(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.linalg.norm(a)
    if denom == 0.0:
        return float(np.linalg.norm(a - b))
    return float(np.linalg.norm(a - b) / denom)


def line_inventory(z: np.ndarray, density: np.ndarray, radius: float = 0.05) -> np.ndarray:
    area = np.pi * radius * radius
    return np.asarray([np.trapz(density[:, ii], z) * area for ii in range(density.shape[1])])


def make_profile_plot(data: dict[str, dict[str, np.ndarray]], out_dir: Path) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharex=True)
    fields = [
        ("n_m", "density [m^-3]"),
        ("Tpar_m", "T_parallel [eV]"),
        ("Tper_m", "T_perp [eV]"),
        ("u_m", "u_parallel [m/s]"),
    ]
    for ax, (field, ylabel) in zip(axes.flat, fields):
        for case_key, meta in CASES.items():
            case = data[case_key]
            if field not in case:
                continue
            ax.plot(case["z"], case[field][:, -1], linewidth=2.0, label=meta["label"])
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.25)
        ax.legend()
    axes[1, 0].set_xlabel("z [m]")
    axes[1, 1].set_xlabel("z [m]")
    fig.suptitle("Final 1 ms hybrid source profiles: main vs PICOS_ECH")
    fig.tight_layout()
    out_path = out_dir / "main_vs_ech_hybrid_final_profiles.png"
    fig.savefig(out_path, dpi=220)
    plt.close(fig)
    return out_path


def make_history_plot(data: dict[str, dict[str, np.ndarray]], out_dir: Path) -> Path:
    fig, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
    for case_key, meta in CASES.items():
        case = data[case_key]
        if "n_m" in case:
            axes[0].plot(case["time"], line_inventory(case["z"], case["n_m"]), linewidth=2.0, label=meta["label"])
        if "Tpar_m" in case and "n_m" in case:
            weights = np.maximum(case["n_m"], 0.0)
            avg_tpar = np.sum(case["Tpar_m"] * weights, axis=0) / np.maximum(np.sum(weights, axis=0), 1.0e-300)
            axes[1].plot(case["time"], avg_tpar, linewidth=2.0, label=meta["label"])
    axes[0].set_ylabel("line inventory [particles]")
    axes[1].set_ylabel("density-weighted T_parallel [eV]")
    axes[1].set_xlabel("time [s]")
    for ax in axes:
        ax.grid(True, alpha=0.25)
        ax.legend()
    fig.suptitle("1 ms hybrid source histories: main vs PICOS_ECH")
    fig.tight_layout()
    out_path = out_dir / "main_vs_ech_hybrid_histories.png"
    fig.savefig(out_path, dpi=220)
    plt.close(fig)
    return out_path


def write_summary(data: dict[str, dict[str, np.ndarray]], out_dir: Path) -> Path:
    rows: list[dict[str, float | str]] = []
    main = data["main"]
    ech = data["ech"]
    for field in ("n_m", "Tpar_m", "Tper_m", "u_m", "BX_m", "EX_m"):
        if field not in main or field not in ech:
            continue
        rows.append(
            {
                "field": field,
                "main_final_min": float(np.nanmin(main[field][:, -1])),
                "main_final_max": float(np.nanmax(main[field][:, -1])),
                "ech_final_min": float(np.nanmin(ech[field][:, -1])),
                "ech_final_max": float(np.nanmax(ech[field][:, -1])),
                "relative_l2_final_ech_minus_main": rel_l2(main[field][:, -1], ech[field][:, -1]),
            }
        )

    out_path = out_dir / "main_vs_ech_hybrid_summary.csv"
    with out_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()) if rows else ["field"])
        writer.writeheader()
        writer.writerows(rows)
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args()

    run_root = args.run_root.resolve()
    out_dir = args.out_dir.resolve() if args.out_dir else run_root / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    data = {key: read_case(hdf_dir(run_root, case)) for key, case in CASES.items()}
    trim_common(data)
    summary = write_summary(data, out_dir)
    profiles = make_profile_plot(data, out_dir)
    histories = make_history_plot(data, out_dir)

    print(f"wrote {summary}")
    print(f"wrote {profiles}")
    print(f"wrote {histories}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
