#!/usr/bin/env python3
"""Lightweight algorithm checks for the kinetic-electron PICOS prototype.

These checks mirror the formulas used by the C++ code and are meant to catch
regressions in the full-orbit Boris update and boundary-aware moment deposition.
They are intentionally dependency-light so they can run on the Mac without
extra Python packages.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass

import numpy as np


@dataclass
class Mesh:
    x_min: float
    x_max: float
    nx: int

    @property
    def dx(self) -> float:
        return (self.x_max - self.x_min) / self.nx

    @property
    def nodes(self) -> np.ndarray:
        return self.x_min + (np.arange(self.nx) + 0.5) * self.dx


def tsc_weights(mesh: Mesh, x: float) -> tuple[int, float, float, float]:
    m = min(max(int(round(0.5 + (x - mesh.x_min) / mesh.dx)) - 1, 0), mesh.nx - 1)
    delta = mesh.nodes[m] - x
    plus = 1.5 + ((delta - mesh.dx) / mesh.dx)
    minus = 1.5 - ((delta + mesh.dx) / mesh.dx)
    wl = 0.5 * plus * plus
    wc = 0.75 - delta * delta / (mesh.dx * mesh.dx)
    wr = 0.5 * minus * minus
    return m, wl, wc, wr


def deposit_particles(mesh: Mesh, x: np.ndarray, periodic: bool) -> np.ndarray:
    scratch = np.zeros(mesh.nx + 4)
    for xp in x:
        m, wl, wc, wr = tsc_weights(mesh, float(xp))
        ix = m + 2
        scratch[ix - 1] += wl
        scratch[ix] += wc
        scratch[ix + 1] += wr

    if periodic:
        scratch[mesh.nx + 1] += scratch[1]
        scratch[2] += scratch[mesh.nx + 2]
    else:
        scratch[2] += scratch[1]
        scratch[mesh.nx + 1] += scratch[mesh.nx + 2]
    scratch[1] = 0.0
    scratch[mesh.nx + 2] = 0.0
    return scratch[2 : mesh.nx + 2]


def boris_step(v: np.ndarray, q_over_m: float, e: np.ndarray, b: np.ndarray, dt: float) -> np.ndarray:
    h = 0.5 * q_over_m * dt
    v_minus = v + h * e
    t = h * b
    s = 2.0 * t / (1.0 + float(np.dot(t, t)))
    v_prime = v_minus + np.cross(v_minus, t)
    v_plus = v_minus + np.cross(v_prime, s)
    return v_plus + h * e


def check_deposition(rng: np.random.Generator, n_particles: int) -> dict[str, float]:
    mesh = Mesh(-5.0, 5.0, 200)
    x = rng.uniform(mesh.x_min, mesh.x_max, size=n_particles)
    # Force shape functions to straddle both boundaries.
    x[:4] = [mesh.x_min, mesh.x_min + 0.1 * mesh.dx, mesh.x_max - 0.1 * mesh.dx, mesh.x_max]

    max_weight_sum_error = 0.0
    min_weight = 1.0
    for xp in x:
        _, wl, wc, wr = tsc_weights(mesh, float(xp))
        max_weight_sum_error = max(max_weight_sum_error, abs((wl + wc + wr) - 1.0))
        min_weight = min(min_weight, wl, wc, wr)

    periodic_total = float(np.sum(deposit_particles(mesh, x, periodic=True)))
    wall_total = float(np.sum(deposit_particles(mesh, x, periodic=False)))
    return {
        "max_weight_sum_error": max_weight_sum_error,
        "min_weight": min_weight,
        "periodic_deposited_total": periodic_total,
        "wall_deposited_total": wall_total,
        "expected_total": float(n_particles),
        "periodic_total_error": abs(periodic_total - n_particles),
        "wall_total_error": abs(wall_total - n_particles),
    }


def check_boris_energy(n_steps: int) -> dict[str, float]:
    q_over_m = -1.0
    e = np.array([0.0, 0.0, 0.0])
    b = np.array([1.0, 0.0, 0.0])
    dt = 0.05
    v = np.array([0.2, 0.7, -0.3])
    initial_energy = 0.5 * float(np.dot(v, v))
    max_relative_error = 0.0
    for _ in range(n_steps):
        v = boris_step(v, q_over_m, e, b, dt)
        energy = 0.5 * float(np.dot(v, v))
        max_relative_error = max(max_relative_error, abs(energy - initial_energy) / initial_energy)
    return {
        "initial_energy": initial_energy,
        "final_energy": 0.5 * float(np.dot(v, v)),
        "max_relative_energy_error": max_relative_error,
        "n_steps": float(n_steps),
    }


def relativistic_speed_from_energy(kinetic_energy: float, mass: float = 1.0, c: float = 1.0) -> float:
    gamma = 1.0 + kinetic_energy / (mass * c * c)
    beta2 = max(0.0, min(1.0 - 1.0e-12, 1.0 - 1.0 / (gamma * gamma)))
    return c * math.sqrt(beta2)


def relativistic_energy_from_speed(speed: float, mass: float = 1.0, c: float = 1.0) -> float:
    beta2 = max(0.0, min(1.0 - 1.0e-12, (speed / c) ** 2))
    gamma = 1.0 / math.sqrt(1.0 - beta2)
    return (gamma - 1.0) * mass * c * c


def check_relativistic_energy_roundtrip() -> dict[str, float]:
    energies = np.array([1.0e-8, 1.0e-4, 1.0e-2, 0.1, 1.0, 10.0, 100.0])
    max_relative_error = 0.0
    max_speed = 0.0
    for energy in energies:
        speed = relativistic_speed_from_energy(float(energy))
        recovered = relativistic_energy_from_speed(speed)
        max_relative_error = max(max_relative_error, abs(recovered - float(energy)) / float(energy))
        max_speed = max(max_speed, float(speed))

    return {
        "max_relative_roundtrip_error": float(max_relative_error),
        "max_speed_over_c": float(max_speed),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--particles", type=int, default=100000)
    parser.add_argument("--boris-steps", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    deposition = check_deposition(rng, args.particles)
    boris = check_boris_energy(args.boris_steps)
    relativistic = check_relativistic_energy_roundtrip()

    ok = bool(
        deposition["max_weight_sum_error"] < 1.0e-12
        and deposition["periodic_total_error"] < 1.0e-9
        and deposition["wall_total_error"] < 1.0e-9
        and deposition["min_weight"] >= -1.0e-12
        and boris["max_relative_energy_error"] < 1.0e-12
        and relativistic["max_relative_roundtrip_error"] < 1.0e-7
        and relativistic["max_speed_over_c"] < 1.0
    )

    print(json.dumps({"ok": ok, "deposition": deposition, "boris": boris, "relativistic": relativistic}, indent=2, sort_keys=True))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
