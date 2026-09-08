#!/usr/bin/env python3
"""Reusable HDF5 reader for PICOS output directories.

The older ProtoLite analysis scripts hard-coded one HDF5 directory and executed
the full plotting workflow at import time. This module keeps only the data
loading pieces and works with either a run directory, a direct HDF5 directory,
or a run tag under ``picosFILES/outputFiles``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import h5py
import numpy as np


E_CHARGE = 1.602176634e-19
M_E = 9.1093837015e-31
C_LIGHT = 299792458.0


def resolve_hdf5_dir(path_or_tag: str | Path, repo_root: Path | None = None) -> Path:
    """Return the HDF5 directory for a PICOS output path or run tag."""

    path = Path(path_or_tag).expanduser()
    candidates = []
    if path.is_absolute():
        candidates.extend([path, path / "HDF5"])
    else:
        root = Path.cwd() if repo_root is None else repo_root
        candidates.extend(
            [
                root / path,
                root / path / "HDF5",
                root / "picosFILES" / "outputFiles" / path / "HDF5",
            ]
        )

    for candidate in candidates:
        if (candidate / "main.h5").is_file():
            return candidate.resolve()
    raise FileNotFoundError(f"Could not find main.h5 for {path_or_tag}")


def _rank_index(path: Path) -> int:
    stem = path.stem
    try:
        return int(stem.rsplit("_", 1)[1])
    except (IndexError, ValueError):
        return 0


def rank_files(hdf5_dir: Path, prefix: str) -> list[Path]:
    files = sorted(hdf5_dir.glob(f"{prefix}_*.h5"), key=_rank_index)
    if not files:
        raise FileNotFoundError(f"No {prefix}_*.h5 files in {hdf5_dir}")
    return files


def numeric_steps(hdf5_dir: Path) -> list[str]:
    """Return ordered output-step group names."""

    for pattern in ("FIELDS_FILE", "PARTICLES_FILE"):
        for path in hdf5_dir.glob(f"{pattern}_*.h5"):
            with h5py.File(path, "r") as handle:
                steps = [key for key in handle.keys() if key.isdigit()]
            if steps:
                return sorted(steps, key=lambda item: int(item))
    raise ValueError(f"No numeric output-step groups found in {hdf5_dir}")


def _scalar(value: Any) -> Any:
    arr = np.asarray(value)
    if arr.shape == ():
        out = arr.item()
    elif arr.size == 1:
        out = arr.reshape(-1)[0].item()
    else:
        return arr
    if isinstance(out, bytes):
        return out.decode()
    return out


def _group_to_dict(group: h5py.Group) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, obj in group.items():
        if isinstance(obj, h5py.Dataset):
            out[key] = _scalar(obj[()])
        elif isinstance(obj, h5py.Group):
            out[key] = _group_to_dict(obj)
    return out


def read_main(hdf5_dir: Path) -> dict[str, Any]:
    with h5py.File(hdf5_dir / "main.h5", "r") as handle:
        return _group_to_dict(handle)


def dataset_exists(file_path: Path, dataset: str) -> bool:
    with h5py.File(file_path, "r") as handle:
        return dataset in handle


def read_dataset(file_path: Path, dataset: str) -> np.ndarray:
    with h5py.File(file_path, "r") as handle:
        return np.asarray(handle[dataset][()])


def available_field_variables(hdf5_dir: Path, step: str | None = None) -> list[str]:
    steps = numeric_steps(hdf5_dir)
    selected = steps[0] if step is None else step
    with h5py.File(rank_files(hdf5_dir, "FIELDS_FILE")[0], "r") as handle:
        root = handle[f"/{selected}/fields"]
        return sorted(root.keys())


def read_field(hdf5_dir: Path, variable: str, step: str) -> np.ndarray:
    parts = []
    dataset = f"/{step}/fields/{variable}/x"
    for path in rank_files(hdf5_dir, "FIELDS_FILE"):
        if dataset_exists(path, dataset):
            parts.append(np.asarray(read_dataset(path, dataset), dtype=float).ravel())
    if not parts:
        raise KeyError(f"Field variable {variable} was not found at step {step}")
    return np.concatenate(parts)


def read_field_series(
    hdf5_dir: Path, variables: Iterable[str] | None = None, steps: Iterable[str] | None = None
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    selected_steps = list(numeric_steps(hdf5_dir) if steps is None else steps)
    available = set(available_field_variables(hdf5_dir, selected_steps[0]))
    selected_vars = sorted(available) if variables is None else [v for v in variables if v in available]

    times = []
    first_field = rank_files(hdf5_dir, "FIELDS_FILE")[0]
    for step in selected_steps:
        times.append(float(_scalar(read_dataset(first_field, f"/{step}/time"))))

    data: dict[str, np.ndarray] = {}
    for variable in selected_vars:
        data[variable] = np.column_stack([read_field(hdf5_dir, variable, step) for step in selected_steps])
    return np.asarray(times, dtype=float), data


def species_names(main: dict[str, Any]) -> list[str]:
    ions = main.get("ions", {})
    n_species = int(ions.get("numberOfParticleSpecies", 0))
    return [f"spp_{idx}" for idx in range(1, n_species + 1)]


def species_metadata(main: dict[str, Any], species: str) -> dict[str, float]:
    ions = main.get("ions", {})
    block = ions.get(species, {})
    return {
        "M": float(block.get("M", np.nan)),
        "Q": float(block.get("Q", np.nan)),
        "Z": float(block.get("Z", np.nan)),
        "NCP": float(block.get("NCP", np.nan)),
        "NSP": float(block.get("NSP", np.nan)),
        "NSP_OUT": float(block.get("NSP_OUT", np.nan)),
        "densityFraction": float(block.get("densityFraction", np.nan)),
    }


def _mesh_dataset_name(variable: str) -> str:
    if variable in {"u_m", "ux_m"}:
        return "u_m/x"
    return variable


def read_species_mesh(hdf5_dir: Path, species: str, variable: str, step: str) -> np.ndarray:
    dataset_name = _mesh_dataset_name(variable)
    dataset = f"/{step}/ions/{species}/{dataset_name}"
    for path in rank_files(hdf5_dir, "PARTICLES_FILE"):
        if dataset_exists(path, dataset):
            return np.asarray(read_dataset(path, dataset), dtype=float).ravel()
    raise KeyError(f"Mesh variable {species}/{variable} was not found at step {step}")


def read_species_mesh_series(
    hdf5_dir: Path,
    species: str,
    variables: Iterable[str],
    steps: Iterable[str] | None = None,
) -> dict[str, np.ndarray]:
    selected_steps = list(numeric_steps(hdf5_dir) if steps is None else steps)
    out: dict[str, np.ndarray] = {}
    for variable in variables:
        columns = []
        for step in selected_steps:
            try:
                columns.append(read_species_mesh(hdf5_dir, species, variable, step))
            except KeyError:
                columns = []
                break
        if columns:
            out[variable] = np.column_stack(columns)
    return out


def velocity_components(velocity: np.ndarray) -> np.ndarray:
    """Return velocity as [component, particle]."""

    velocity = np.asarray(velocity, dtype=float)
    if velocity.ndim != 2:
        raise ValueError(f"Expected 2D V_p array, got shape {velocity.shape}")
    if velocity.shape[0] in (2, 3):
        return velocity
    if velocity.shape[1] in (2, 3):
        return velocity.T
    raise ValueError(f"Cannot infer velocity component axis from shape {velocity.shape}")


def read_species_particle(hdf5_dir: Path, species: str, variable: str, step: str, stride: int = 1) -> np.ndarray:
    dataset = f"/{step}/ions/{species}/{variable}"
    arrays = []
    for path in rank_files(hdf5_dir, "PARTICLES_FILE"):
        if not dataset_exists(path, dataset):
            continue
        values = np.asarray(read_dataset(path, dataset), dtype=float)
        if variable == "V_p":
            values = velocity_components(values)
            values = values[:, ::stride]
            arrays.append(values)
        else:
            arrays.append(values.ravel()[::stride])
    if not arrays:
        raise KeyError(f"Particle variable {species}/{variable} was not found at step {step}")
    if variable == "V_p":
        return np.concatenate(arrays, axis=1)
    return np.concatenate(arrays)


def read_species_particle_state(
    hdf5_dir: Path,
    species: str,
    step: str,
    variables: Iterable[str] = ("X_p", "V_p", "a_p", "BX_p", "EX_p"),
    stride: int = 1,
    max_particles: int | None = None,
) -> dict[str, np.ndarray]:
    state: dict[str, np.ndarray] = {}
    for variable in variables:
        try:
            values = read_species_particle(hdf5_dir, species, variable, step, stride=stride)
        except KeyError:
            continue
        if max_particles is not None:
            if variable == "V_p":
                values = values[:, :max_particles]
            else:
                values = values[:max_particles]
        state[variable] = values
    return state


def kinetic_energy_eV(velocity: np.ndarray, mass_kg: float, relativistic: bool = False) -> np.ndarray:
    components = velocity_components(velocity)
    speed2 = np.sum(components * components, axis=0)
    if not relativistic:
        return 0.5 * mass_kg * speed2 / E_CHARGE
    beta2 = np.clip(speed2 / (C_LIGHT * C_LIGHT), 0.0, 1.0 - 1.0e-15)
    gamma = 1.0 / np.sqrt(1.0 - beta2)
    return (gamma - 1.0) * mass_kg * C_LIGHT * C_LIGHT / E_CHARGE


def weighted_mean(values: np.ndarray, weights: np.ndarray | None = None) -> float:
    values = np.asarray(values, dtype=float)
    mask = np.isfinite(values)
    if weights is None:
        selected = values[mask]
        return float(np.mean(selected)) if selected.size else float("nan")
    weights = np.asarray(weights, dtype=float)
    mask &= np.isfinite(weights) & (weights > 0.0)
    if not np.any(mask):
        return float("nan")
    return float(np.average(values[mask], weights=weights[mask]))


def weighted_percentile(values: np.ndarray, percentile: float, weights: np.ndarray | None = None) -> float:
    values = np.asarray(values, dtype=float)
    mask = np.isfinite(values)
    if weights is None:
        selected = values[mask]
        return float(np.percentile(selected, percentile)) if selected.size else float("nan")
    weights = np.asarray(weights, dtype=float)
    mask &= np.isfinite(weights) & (weights > 0.0)
    if not np.any(mask):
        return float("nan")
    x = values[mask]
    w = weights[mask]
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    cumulative = np.cumsum(w)
    target = percentile / 100.0 * cumulative[-1]
    return float(np.interp(target, cumulative, x))


def rolling_mean(values: np.ndarray, window: int, axis: int = 0) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if window <= 1:
        return values
    kernel = np.ones(window, dtype=float) / float(window)
    return np.apply_along_axis(lambda line: np.convolve(line, kernel, mode="same"), axis, values)


def central_diff(values: np.ndarray, spacing: float, axis: int = 0) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return np.gradient(values, spacing, axis=axis, edge_order=1)


@dataclass(frozen=True)
class PicosRun:
    hdf5_dir: Path
    main: dict[str, Any]
    steps: list[str]

    @classmethod
    def open(cls, path_or_tag: str | Path, repo_root: Path | None = None) -> "PicosRun":
        hdf5_dir = resolve_hdf5_dir(path_or_tag, repo_root=repo_root)
        return cls(hdf5_dir=hdf5_dir, main=read_main(hdf5_dir), steps=numeric_steps(hdf5_dir))

    @property
    def x_m(self) -> np.ndarray:
        geometry = self.main.get("geometry", {})
        if "x_m" in geometry:
            return np.asarray(geometry["x_m"], dtype=float)
        if "xAxis" in geometry:
            return np.asarray(geometry["xAxis"], dtype=float)
        raise KeyError("main.h5 does not contain geometry/x_m or geometry/xAxis")

    @property
    def species(self) -> list[str]:
        return species_names(self.main)

