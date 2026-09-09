#!/usr/bin/env python3
"""Analytic checks for the PICOS++ Coulomb collision operator.

The tests here are intentionally small and deterministic.  They validate the
same scalar formulas used by `collisionOperator.cpp` and the first moments of
the Monte Carlo pitch/energy updates against constant-coefficient analytic
solutions:

    <xi(t)> = xi0 exp(-nu_D t)
    <E(t)>  = Eeq + (E0 - Eeq) exp(-2 nu_E t)

The constant-coefficient checks isolate the stochastic update math from the
velocity-dependent Coulomb coefficients.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


F_E = 1.602176e-19
F_ME = 9.109382e-31
F_MP = 1.672621e-27
F_EPSILON = 8.854e-12


def phi(x: np.ndarray) -> np.ndarray:
    return np.vectorize(math.erf)(x)


def phip(x: np.ndarray) -> np.ndarray:
    return 2.0 / math.sqrt(math.pi) * np.exp(-x * x)


def gb(x: np.ndarray) -> np.ndarray:
    out = np.empty_like(x)
    small = x < 0.01
    out[small] = (2.0 / math.sqrt(math.pi) / 3.0) * x[small]
    out[~small] = (phi(x[~small]) - x[~small] * phip(x[~small])) / (2.0 * x[~small] ** 2)
    return out


def e_nue_derivative(x: np.ndarray) -> np.ndarray:
    out = np.zeros_like(x)
    mask = x >= 0.01
    denominator = phi(x[mask]) - x[mask] * phip(x[mask])
    phi2p = -4.0 * x[mask] / math.sqrt(math.pi) * np.exp(-x[mask] * x[mask])
    out[mask] = 0.5 * (
        (3.0 * x[mask] * phip(x[mask]) - 3.0 * phi(x[mask]) - x[mask] ** 2 * phi2p)
        / denominator
    )
    return out


def log_lambda(nb_m3: float, tb_ev: float) -> float:
    return 30.0 - 0.5 * math.log(nb_m3 / (tb_ev ** 1.5))


def nu_ab0(nb_m3: float, tb_ev: float, mb_kg: float, za: float, zb: float, ma_kg: float) -> float:
    wtb = math.sqrt(2.0 * F_E * tb_ev / mb_kg)
    return (
        nb_m3
        * F_E**4
        * (za * zb) ** 2
        * log_lambda(nb_m3, tb_ev)
        / (2.0 * math.pi * ma_kg * ma_kg * F_EPSILON * F_EPSILON * wtb**3)
    )


def nu_d(x: np.ndarray, nu0: float) -> np.ndarray:
    return nu0 * (phi(x) - gb(x)) / (x**3)


def nu_e_boozer_kim(x: np.ndarray, nu0: float, ma_kg: float, mb_kg: float) -> np.ndarray:
    return nu0 * (2.0 * ma_kg / mb_kg) * gb(x) / x / (1.0 + mb_kg / ma_kg)


def reflect_pitch(xi: np.ndarray) -> np.ndarray:
    over = xi * xi > 1.0
    xi[over] = np.sign(xi[over]) - np.fmod(xi[over], np.sign(xi[over]))
    return xi


def validate_rates(out_dir: Path) -> dict[str, float]:
    x = np.geomspace(1.0e-3, 20.0, 500)
    ma = F_ME
    mb = F_ME
    nu0 = nu_ab0(1.0e18, 15.0, mb, -1.0, -1.0, ma)

    picos_nud = nu_d(x, nu0)
    analytic_nud = nu0 * (phi(x) - gb(x)) / (x**3)
    picos_nue = nu_e_boozer_kim(x, nu0, ma, mb)
    analytic_nue = nu0 * (2.0 * ma / mb) * gb(x) / x / (1.0 + mb / ma)

    rel_nud = np.max(np.abs(picos_nud - analytic_nud) / np.maximum(np.abs(analytic_nud), 1.0e-300))
    rel_nue = np.max(np.abs(picos_nue - analytic_nue) / np.maximum(np.abs(analytic_nue), 1.0e-300))

    with (out_dir / "collision_rate_formula_check.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x", "nu_D_s-1", "nu_E_boozer_kim_s-1", "E_nuE_d_nuE_dE"])
        for row in zip(x, picos_nud, picos_nue, e_nue_derivative(x)):
            writer.writerow([f"{value:.16e}" for value in row])

    return {
        "rate_nu_D_max_relative_error": float(rel_nud),
        "rate_nu_E_max_relative_error": float(rel_nue),
    }


def validate_pitch_moment(out_dir: Path, rng: np.random.Generator) -> dict[str, float]:
    particles = 120_000
    nu = 2.0e4
    dt = 1.0e-7
    steps = 700
    xi0 = 0.75
    xi = np.full(particles, xi0)

    times = np.empty(steps + 1)
    means = np.empty(steps + 1)
    analytic = np.empty(steps + 1)
    nudt = nu * dt
    for step in range(steps + 1):
        t = step * dt
        times[step] = t
        means[step] = float(np.mean(xi))
        analytic[step] = xi0 * math.exp(-nu * t)
        if step == steps:
            break
        signs = rng.choice(np.array([-1.0, 1.0]), size=particles)
        xi = xi * (1.0 - nudt) + signs * np.sqrt(np.maximum(0.0, (1.0 - xi * xi) * nudt))
        xi = reflect_pitch(xi)

    err = np.abs(means - analytic)
    max_abs = float(np.max(err))
    final_abs = float(err[-1])
    with (out_dir / "pitch_first_moment.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_s", "mean_xi_mc", "mean_xi_analytic"])
        for row in zip(times, means, analytic):
            writer.writerow([f"{value:.16e}" for value in row])

    return {
        "pitch_mean_max_abs_error": max_abs,
        "pitch_mean_final_abs_error": final_abs,
    }


def validate_energy_moment(out_dir: Path, rng: np.random.Generator) -> dict[str, float]:
    particles = 120_000
    nu = 1.0e4
    dt = 1.0e-7
    steps = 900
    tb_ev = 15.0
    e_eq = 1.5 * tb_ev
    e0 = 120.0
    energy = np.full(particles, e0)

    times = np.empty(steps + 1)
    means = np.empty(steps + 1)
    analytic = np.empty(steps + 1)
    nudt = nu * dt
    for step in range(steps + 1):
        t = step * dt
        times[step] = t
        means[step] = float(np.mean(energy))
        analytic[step] = e_eq + (e0 - e_eq) * math.exp(-2.0 * nu * t)
        if step == steps:
            break
        signs = rng.choice(np.array([-2.0, 2.0]), size=particles)
        energy = (
            energy * (1.0 - 2.0 * nudt)
            + 2.0 * nudt * e_eq
            + signs * np.sqrt(np.maximum(0.0, tb_ev * nudt * energy))
        )
        energy = np.maximum(energy, 0.0)

    err = np.abs(means - analytic)
    max_rel = float(np.max(err / np.maximum(np.abs(analytic), 1.0e-12)))
    final_rel = float(err[-1] / max(abs(analytic[-1]), 1.0e-12))
    with (out_dir / "energy_first_moment.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_s", "mean_energy_mc_ev", "mean_energy_analytic_ev"])
        for row in zip(times, means, analytic):
            writer.writerow([f"{value:.16e}" for value in row])

    return {
        "energy_mean_max_relative_error": max_rel,
        "energy_mean_final_relative_error": final_rel,
    }


def write_plots(out_dir: Path) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from paper_plot_style import apply_paper_figure_style

        apply_paper_figure_style()
    except Exception:
        return

    pitch = np.genfromtxt(out_dir / "pitch_first_moment.csv", delimiter=",", names=True)
    energy = np.genfromtxt(out_dir / "energy_first_moment.csv", delimiter=",", names=True)
    rates = np.genfromtxt(out_dir / "collision_rate_formula_check.csv", delimiter=",", names=True)

    fig, axes = plt.subplots(1, 3, figsize=(21, 6), constrained_layout=True)
    axes[0].plot(pitch["time_s"], pitch["mean_xi_mc"], label="Monte Carlo")
    axes[0].plot(pitch["time_s"], pitch["mean_xi_analytic"], "--", label="analytic")
    axes[0].set_xlabel("time [s]")
    axes[0].set_ylabel(r"$\langle \xi \rangle$")
    axes[0].set_title("Pitch Moment")
    axes[0].legend()

    axes[1].plot(energy["time_s"], energy["mean_energy_mc_ev"], label="Monte Carlo")
    axes[1].plot(energy["time_s"], energy["mean_energy_analytic_ev"], "--", label="analytic")
    axes[1].set_xlabel("time [s]")
    axes[1].set_ylabel(r"$\langle E \rangle$ [eV]")
    axes[1].set_title("Energy Moment")
    axes[1].legend()

    axes[2].loglog(rates["x"], rates["nu_D_s1"], label="nu_D")
    axes[2].loglog(rates["x"], rates["nu_E_boozer_kim_s1"], label="nu_E Boozer-Kim")
    axes[2].set_xlabel("x = v/v_Tb")
    axes[2].set_ylabel(r"rate [s$^{-1}$]")
    axes[2].set_title("Coulomb Rates")
    axes[2].legend()

    fig.savefig(out_dir / "collision_operator_analytic_validation.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("validation/collision_operator_analytic"))
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--pitch-max-abs-tol", type=float, default=8.0e-3)
    parser.add_argument("--energy-max-rel-tol", type=float, default=2.5e-2)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    summary = {}
    summary.update(validate_rates(args.output_dir))
    summary.update(validate_pitch_moment(args.output_dir, rng))
    summary.update(validate_energy_moment(args.output_dir, rng))
    summary["passed"] = (
        summary["rate_nu_D_max_relative_error"] < 1.0e-14
        and summary["rate_nu_E_max_relative_error"] < 1.0e-14
        and summary["pitch_mean_max_abs_error"] < args.pitch_max_abs_tol
        and summary["energy_mean_max_relative_error"] < args.energy_max_rel_tol
    )

    write_plots(args.output_dir)
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
