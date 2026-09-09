#!/usr/bin/env python3
"""Analytic checks for PICOS++ electrostatic fields and particle sources.

The tests are deterministic and lightweight.  They mirror the one-dimensional
finite-difference/FFT Poisson updates and the particle source sampling formulas
used in the PICOS++ kinetic-electron branch.

Validated pieces:

* Dirichlet charge-density Poisson solve using a manufactured sine potential.
* Periodic charge-density Poisson solve using a manufactured cosine potential.
* Reformulated-Poisson stress-moment E-field update using a manufactured stress.
* Coupled pair-source Gaussian spatial sampling against a truncated normal.
* Coupled pair-source profile sampling against the expected discrete weights.
* Pair-source ion/electron weights against charge-neutral injection.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


EPS0 = 8.8541878128e-12
M_E_OVER_M_D = 9.1093837015e-31 / (2.0 * 1.66053906660e-27)


def normal_pdf(x: np.ndarray | float) -> np.ndarray | float:
    return np.exp(-0.5 * np.asarray(x) ** 2) / math.sqrt(2.0 * math.pi)


def normal_cdf(x: np.ndarray | float) -> np.ndarray | float:
    vec_erf = np.vectorize(math.erf)
    return 0.5 * (1.0 + vec_erf(np.asarray(x) / math.sqrt(2.0)))


def thomas_solve(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    lower = lower.copy()
    diag = diag.copy()
    upper = upper.copy()
    rhs = rhs.copy()
    n = rhs.size
    for i in range(1, n):
        factor = lower[i] / diag[i - 1]
        diag[i] -= factor * upper[i - 1]
        rhs[i] -= factor * rhs[i - 1]
    out = np.zeros(n)
    out[-1] = rhs[-1] / diag[-1]
    for i in range(n - 2, -1, -1):
        out[i] = (rhs[i] - upper[i] * out[i + 1]) / diag[i]
    return out


def solve_dirichlet_poisson(rho: np.ndarray, dx: float, phi_left: float, phi_right: float, eps: float = EPS0) -> np.ndarray:
    n = rho.size
    lower = np.ones(n)
    diag = -2.0 * np.ones(n)
    upper = np.ones(n)
    lower[0] = 0.0
    upper[-1] = 0.0
    rhs = -rho * dx * dx / eps
    rhs[0] -= phi_left
    rhs[-1] -= phi_right
    phi = np.zeros(n + 2)
    phi[0] = phi_left
    phi[-1] = phi_right
    phi[1:-1] = thomas_solve(lower, diag, upper, rhs)
    return phi


def solve_periodic_poisson(rho: np.ndarray, dx: float, eps: float = EPS0) -> np.ndarray:
    n = rho.size
    rhs = -rho * dx * dx / eps
    rhs -= np.mean(rhs)
    rhs_hat = np.fft.fft(rhs)
    phi_hat = np.zeros(n, dtype=complex)
    for kk in range(1, n):
        eig = -4.0 * math.sin(math.pi * kk / n) ** 2
        phi_hat[kk] = rhs_hat[kk] / eig
    solution = np.real(np.fft.ifft(phi_hat))
    phi = np.zeros(n + 2)
    phi[1:-1] = solution
    phi[0] = phi[-2]
    phi[-1] = phi[1]
    return phi


def electric_from_phi(phi: np.ndarray, dx: float, periodic: bool) -> np.ndarray:
    n = phi.size - 2
    e = np.zeros(n + 2)
    e[1:-1] = -(phi[2:] - phi[:-2]) / (2.0 * dx)
    if periodic:
        e[0] = e[-2]
        e[-1] = e[1]
    else:
        e[0] = e[1]
        e[-1] = e[-2]
    return e


def l2_error(numerical: np.ndarray, exact: np.ndarray) -> float:
    return float(np.sqrt(np.mean((numerical - exact) ** 2)))


def convergence_order(errors: list[float]) -> float:
    if len(errors) < 2 or errors[-1] <= 0.0 or errors[-2] <= 0.0:
        return math.nan
    return float(math.log(errors[-2] / errors[-1]) / math.log(2.0))


def validate_dirichlet_poisson(out_dir: Path, ns: list[int]) -> tuple[list[dict[str, float]], dict[str, np.ndarray]]:
    rows: list[dict[str, float]] = []
    last_plot: dict[str, np.ndarray] = {}
    length = 1.0
    amplitude = 3.0
    for n in ns:
        dx = length / (n + 1)
        x = dx * np.arange(1, n + 1)
        k = math.pi / length
        phi_exact = amplitude * np.sin(k * x)
        rho = EPS0 * amplitude * k * k * np.sin(k * x)
        phi = solve_dirichlet_poisson(rho, dx, 0.0, 0.0)
        e = electric_from_phi(phi, dx, periodic=False)
        e_exact = -amplitude * k * np.cos(k * x)
        rows.append(
            {
                "N": float(n),
                "dx": dx,
                "phi_l2_error": l2_error(phi[1:-1], phi_exact),
                "phi_max_error": float(np.max(np.abs(phi[1:-1] - phi_exact))),
                "E_l2_error": l2_error(e[1:-1], e_exact),
                "E_max_error": float(np.max(np.abs(e[1:-1] - e_exact))),
            }
        )
        if n == ns[-2] if len(ns) > 1 else n == ns[-1]:
            last_plot = {
                "x": x,
                "phi_numeric": phi[1:-1],
                "phi_exact": phi_exact,
                "E_numeric": e[1:-1],
                "E_exact": e_exact,
            }
    return rows, last_plot


def validate_periodic_poisson(out_dir: Path, ns: list[int]) -> tuple[list[dict[str, float]], dict[str, np.ndarray]]:
    rows: list[dict[str, float]] = []
    last_plot: dict[str, np.ndarray] = {}
    length = 1.0
    amplitude = 2.0
    mode = 2
    for n in ns:
        dx = length / n
        x = dx * np.arange(n)
        k = 2.0 * math.pi * mode / length
        phi_exact = amplitude * np.cos(k * x)
        rho = EPS0 * amplitude * k * k * np.cos(k * x)
        phi = solve_periodic_poisson(rho, dx)
        phi_num = phi[1:-1] - np.mean(phi[1:-1])
        e = electric_from_phi(phi, dx, periodic=True)
        e_exact = amplitude * k * np.sin(k * x)
        rows.append(
            {
                "N": float(n),
                "dx": dx,
                "phi_l2_error": l2_error(phi_num, phi_exact),
                "phi_max_error": float(np.max(np.abs(phi_num - phi_exact))),
                "E_l2_error": l2_error(e[1:-1], e_exact),
                "E_max_error": float(np.max(np.abs(e[1:-1] - e_exact))),
            }
        )
        if n == ns[-2] if len(ns) > 1 else n == ns[-1]:
            last_plot = {
                "x": x,
                "phi_numeric": phi_num,
                "phi_exact": phi_exact,
                "E_numeric": e[1:-1],
                "E_exact": e_exact,
            }
    return rows, last_plot


def validate_reformulated_poisson(out_dir: Path) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    n = 512
    length = 1.0
    dx = length / n
    x = dx * np.arange(n)
    k = 2.0 * math.pi / length
    e0 = 4.0
    density_denominator = 1.0 + 1.0 / M_E_OVER_M_D
    e_exact = e0 * np.sin(k * x)

    stress = -density_denominator * e0 * np.cos(k * x) / k
    stress_g = np.zeros(n + 2)
    stress_g[1:-1] = stress
    stress_g[0] = stress_g[-2]
    stress_g[-1] = stress_g[1]
    div_stress = (stress_g[2:] - stress_g[:-2]) / (2.0 * dx)

    e_quasineutral = div_stress / density_denominator
    lam = 0.02
    dt = 1.0e-5
    alpha = (dt * density_denominator / (lam * lam)) / (1.0 + dt * density_denominator / (lam * lam))
    e_implicit = (dt * div_stress / (lam * lam)) / (1.0 + dt * density_denominator / (lam * lam))
    e_implicit_exact = alpha * e_exact

    row = {
        "N": float(n),
        "dx": dx,
        "quasineutral_E_l2_error": l2_error(e_quasineutral, e_exact),
        "quasineutral_E_max_error": float(np.max(np.abs(e_quasineutral - e_exact))),
        "implicit_E_l2_error": l2_error(e_implicit, e_implicit_exact),
        "implicit_E_max_error": float(np.max(np.abs(e_implicit - e_implicit_exact))),
        "implicit_alpha": float(alpha),
    }
    plot = {
        "x": x,
        "E_exact": e_exact,
        "E_quasineutral": e_quasineutral,
        "E_implicit": e_implicit,
        "E_implicit_exact": e_implicit_exact,
    }
    return row, plot


def sample_gaussian_source(
    rng: np.random.Generator,
    count: int,
    mean_x: float,
    sigma_x: float,
    xmin: float,
    xmax: float,
) -> np.ndarray:
    out = np.empty(count)
    filled = 0
    sigma = sigma_x * math.sqrt(2.0)
    while filled < count:
        batch = max(8192, count - filled)
        u = np.maximum(rng.random(batch), np.finfo(float).tiny)
        theta = 2.0 * math.pi * rng.random(batch)
        samples = mean_x + sigma * np.sqrt(-np.log(u)) * np.cos(theta)
        samples = samples[(samples >= xmin) & (samples <= xmax)]
        take = min(samples.size, count - filled)
        if take:
            out[filled : filled + take] = samples[:take]
            filled += take
    return out


def truncated_normal_moments(mean_x: float, sigma_x: float, xmin: float, xmax: float) -> tuple[float, float]:
    alpha = (xmin - mean_x) / sigma_x
    beta = (xmax - mean_x) / sigma_x
    z_norm = float(normal_cdf(beta) - normal_cdf(alpha))
    mean = mean_x + float((normal_pdf(alpha) - normal_pdf(beta)) / z_norm) * sigma_x
    variance = sigma_x * sigma_x * (
        1.0
        + float((alpha * normal_pdf(alpha) - beta * normal_pdf(beta)) / z_norm)
        - float((normal_pdf(alpha) - normal_pdf(beta)) / z_norm) ** 2
    )
    return mean, math.sqrt(max(variance, 0.0))


def truncated_normal_cdf(x: np.ndarray, mean_x: float, sigma_x: float, xmin: float, xmax: float) -> np.ndarray:
    alpha = (xmin - mean_x) / sigma_x
    beta = (xmax - mean_x) / sigma_x
    norm = normal_cdf(beta) - normal_cdf(alpha)
    val = (normal_cdf((x - mean_x) / sigma_x) - normal_cdf(alpha)) / norm
    return np.clip(val, 0.0, 1.0)


def validate_gaussian_pair_source(out_dir: Path, rng: np.random.Generator) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    count = 200_000
    xmin = -0.05
    xmax = 0.05
    mean_x = 0.012
    sigma_x = 0.018
    samples = sample_gaussian_source(rng, count, mean_x, sigma_x, xmin, xmax)
    analytic_mean, analytic_std = truncated_normal_moments(mean_x, sigma_x, xmin, xmax)
    sorted_samples = np.sort(samples)
    empirical = np.arange(1, count + 1) / count
    analytic = truncated_normal_cdf(sorted_samples, mean_x, sigma_x, xmin, xmax)
    cdf_error = float(np.max(np.abs(empirical - analytic)))
    hist, edges = np.histogram(samples, bins=90, range=(xmin, xmax), density=True)
    centers = 0.5 * (edges[1:] + edges[:-1])
    pdf = normal_pdf((centers - mean_x) / sigma_x) / sigma_x
    pdf /= float(normal_cdf((xmax - mean_x) / sigma_x) - normal_cdf((xmin - mean_x) / sigma_x))
    row = {
        "samples": float(count),
        "sample_mean_m": float(np.mean(samples)),
        "analytic_mean_m": analytic_mean,
        "sample_std_m": float(np.std(samples)),
        "analytic_std_m": analytic_std,
        "max_cdf_error": cdf_error,
    }
    plot = {
        "centers": centers,
        "hist_pdf": hist,
        "analytic_pdf": pdf,
        "samples": samples,
    }
    return row, plot


def validate_profile_pair_source(out_dir: Path, rng: np.random.Generator) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    count = 300_000
    xmin = -0.05
    xmax = 0.05
    n_profile = 80
    x_profile = np.linspace(xmin, xmax, n_profile)
    length = xmax - xmin
    profile = (
        1.0
        + 0.45 * np.cos(2.0 * math.pi * (x_profile - xmin) / length)
        + 0.65 * np.exp(-((x_profile - 0.018) / 0.011) ** 2)
    )
    weights = np.maximum(profile, 0.0)
    expected = weights / np.sum(weights)
    selected = rng.choice(np.arange(n_profile), size=count, p=expected)
    observed = np.bincount(selected, minlength=n_profile) / count
    dx = abs(x_profile[1] - x_profile[0])
    jittered = np.clip(x_profile[selected] + (rng.random(count) - 0.5) * dx, xmin, xmax)
    row = {
        "samples": float(count),
        "profile_cells": float(n_profile),
        "l1_probability_error": float(np.sum(np.abs(observed - expected))),
        "max_probability_error": float(np.max(np.abs(observed - expected))),
        "sample_mean_m": float(np.mean(jittered)),
        "expected_mean_m": float(np.sum(x_profile * expected)),
    }
    plot = {
        "x_profile": x_profile,
        "expected": expected,
        "observed": observed,
        "jittered": jittered,
    }
    return row, plot


def pair_source_weights(
    rate: float,
    dt: float,
    global_pairs: float,
    ion_ncp: float,
    electron_ncp: float,
    ion_z: float,
    electron_z: float,
    max_weight: float,
) -> tuple[float, float]:
    real_ion_pairs_per_step = rate * dt
    ion_weight = min(real_ion_pairs_per_step / (ion_ncp * global_pairs), max_weight)
    electron_weight = min(
        real_ion_pairs_per_step
        * abs(ion_z)
        / (abs(electron_z) * electron_ncp * global_pairs),
        max_weight,
    )
    return ion_weight, electron_weight


def validate_pair_source_charge_neutrality(out_dir: Path) -> list[dict[str, float]]:
    cases = [
        {
            "name": "D_plus_electron",
            "rate": 1.0e19,
            "dt": 2.0e-8,
            "pairs": 5000.0,
            "ion_ncp": 1.0e8,
            "electron_ncp": 1.0e8,
            "ion_z": 1.0,
            "electron_z": -1.0,
            "max_weight": 1000.0,
        },
        {
            "name": "Z2_ion_electron_unequal_NCP",
            "rate": 4.0e18,
            "dt": 1.0e-8,
            "pairs": 4000.0,
            "ion_ncp": 2.0e8,
            "electron_ncp": 5.0e7,
            "ion_z": 2.0,
            "electron_z": -1.0,
            "max_weight": 1000.0,
        },
    ]
    rows: list[dict[str, float]] = []
    for case in cases:
        ion_w, electron_w = pair_source_weights(
            case["rate"],
            case["dt"],
            case["pairs"],
            case["ion_ncp"],
            case["electron_ncp"],
            case["ion_z"],
            case["electron_z"],
            case["max_weight"],
        )
        ion_real = case["pairs"] * case["ion_ncp"] * ion_w
        electron_real = case["pairs"] * case["electron_ncp"] * electron_w
        net_charge_units = case["ion_z"] * ion_real + case["electron_z"] * electron_real
        denom = abs(case["ion_z"] * ion_real) + abs(case["electron_z"] * electron_real)
        rows.append(
            {
                "case": case["name"],
                "ion_weight": ion_w,
                "electron_weight": electron_w,
                "ion_real_particles": ion_real,
                "electron_real_particles": electron_real,
                "net_charge_normalized": float(net_charge_units / denom),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_plot(out_dir: Path, dirichlet: dict[str, np.ndarray], periodic: dict[str, np.ndarray],
               reformulated: dict[str, np.ndarray], gaussian: dict[str, np.ndarray],
               profile: dict[str, np.ndarray], summary: dict[str, float]) -> Path | None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from paper_plot_style import apply_paper_figure_style

        apply_paper_figure_style()
    except Exception:
        return None

    fig, axes = plt.subplots(2, 3, figsize=(21, 12), constrained_layout=True)

    ax = axes[0, 0]
    ax.plot(dirichlet["x"], dirichlet["phi_exact"], "k-", label="analytic")
    ax.plot(dirichlet["x"], dirichlet["phi_numeric"], "r--", label="numeric")
    ax.set_title("Dirichlet Poisson")
    ax.set_xlabel("z [m]")
    ax.set_ylabel(r"$\phi$ [V]")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.plot(periodic["x"], periodic["phi_exact"], "k-", label="analytic")
    ax.plot(periodic["x"], periodic["phi_numeric"], "b--", label="numeric")
    ax.set_title("Periodic Poisson")
    ax.set_xlabel("z [m]")
    ax.set_ylabel(r"$\phi$ [V]")
    ax.legend(frameon=False)

    ax = axes[0, 2]
    ax.plot(reformulated["x"], reformulated["E_exact"], "k-", label="analytic")
    ax.plot(reformulated["x"], reformulated["E_quasineutral"], "g--", label="quasi-neutral")
    ax.plot(reformulated["x"], reformulated["E_implicit"], "m:", label="implicit")
    ax.set_title("Reformulated Poisson")
    ax.set_xlabel("z [m]")
    ax.set_ylabel(r"$E_\parallel$ [arb.]")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.semilogy(
        summary["poisson_dx"],
        summary["dirichlet_phi_l2"],
        "o-",
        label="Dirichlet",
    )
    ax.semilogy(
        summary["periodic_dx"],
        summary["periodic_phi_l2"],
        "s-",
        label="Periodic",
    )
    ax.invert_xaxis()
    ax.set_title("Poisson convergence")
    ax.set_xlabel("grid spacing [m]")
    ax.set_ylabel(r"$L_2(\phi)$")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ax.plot(gaussian["centers"], gaussian["hist_pdf"], color="tab:orange", label="sampled")
    ax.plot(gaussian["centers"], gaussian["analytic_pdf"], "k--", label="truncated normal")
    ax.set_title("Gaussian pair source")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("probability density")
    ax.legend(frameon=False)

    ax = axes[1, 2]
    ax.plot(profile["x_profile"], profile["expected"], "k-", label="expected")
    ax.plot(profile["x_profile"], profile["observed"], "c.", label="sampled")
    ax.set_title("Profile pair source")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("cell probability")
    ax.legend(frameon=False)

    path = out_dir / "field_source_analytic.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def assert_validation(summary: dict[str, float]) -> None:
    checks = [
        ("dirichlet_phi_order", summary["dirichlet_phi_order"], 1.8, None),
        ("periodic_phi_order", summary["periodic_phi_order"], 1.8, None),
        ("reformulated_quasineutral_E_max_error", summary["reformulated_quasineutral_E_max_error"], None, 5.0e-4),
        ("gaussian_source_max_cdf_error", summary["gaussian_source_max_cdf_error"], None, 5.0e-3),
        ("profile_source_l1_probability_error", summary["profile_source_l1_probability_error"], None, 4.0e-2),
        ("pair_source_max_abs_net_charge_normalized", summary["pair_source_max_abs_net_charge_normalized"], None, 1.0e-14),
    ]
    failures = []
    for name, value, lower, upper in checks:
        if lower is not None and value < lower:
            failures.append(f"{name}={value:.6e} < {lower:.6e}")
        if upper is not None and value > upper:
            failures.append(f"{name}={value:.6e} > {upper:.6e}")
    if failures:
        raise SystemExit("validation failed: " + "; ".join(failures))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("validation/field_source_analytic"))
    parser.add_argument("--seed", type=int, default=20260909)
    parser.add_argument("--assert-validation", action="store_true")
    args = parser.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    ns = [32, 64, 128, 256, 512]
    dirichlet_rows, dirichlet_plot = validate_dirichlet_poisson(out_dir, ns)
    periodic_rows, periodic_plot = validate_periodic_poisson(out_dir, ns)
    reformulated_row, reformulated_plot = validate_reformulated_poisson(out_dir)
    gaussian_row, gaussian_plot = validate_gaussian_pair_source(out_dir, rng)
    profile_row, profile_plot = validate_profile_pair_source(out_dir, rng)
    charge_rows = validate_pair_source_charge_neutrality(out_dir)

    write_csv(out_dir / "poisson_dirichlet_convergence.csv", dirichlet_rows)
    write_csv(out_dir / "poisson_periodic_convergence.csv", periodic_rows)
    write_csv(out_dir / "reformulated_poisson_check.csv", [reformulated_row])
    write_csv(out_dir / "pair_source_gaussian_check.csv", [gaussian_row])
    write_csv(out_dir / "pair_source_profile_check.csv", [profile_row])
    write_csv(out_dir / "pair_source_charge_neutrality.csv", charge_rows)

    summary: dict[str, float | str | list[float]] = {
        "seed": float(args.seed),
        "dirichlet_phi_order": convergence_order([row["phi_l2_error"] for row in dirichlet_rows]),
        "dirichlet_E_order": convergence_order([row["E_l2_error"] for row in dirichlet_rows]),
        "dirichlet_phi_l2_finest": dirichlet_rows[-1]["phi_l2_error"],
        "dirichlet_E_l2_finest": dirichlet_rows[-1]["E_l2_error"],
        "periodic_phi_order": convergence_order([row["phi_l2_error"] for row in periodic_rows]),
        "periodic_E_order": convergence_order([row["E_l2_error"] for row in periodic_rows]),
        "periodic_phi_l2_finest": periodic_rows[-1]["phi_l2_error"],
        "periodic_E_l2_finest": periodic_rows[-1]["E_l2_error"],
        "reformulated_quasineutral_E_max_error": reformulated_row["quasineutral_E_max_error"],
        "reformulated_implicit_E_max_error": reformulated_row["implicit_E_max_error"],
        "gaussian_source_mean_abs_error_m": abs(gaussian_row["sample_mean_m"] - gaussian_row["analytic_mean_m"]),
        "gaussian_source_std_abs_error_m": abs(gaussian_row["sample_std_m"] - gaussian_row["analytic_std_m"]),
        "gaussian_source_max_cdf_error": gaussian_row["max_cdf_error"],
        "profile_source_l1_probability_error": profile_row["l1_probability_error"],
        "profile_source_max_probability_error": profile_row["max_probability_error"],
        "pair_source_max_abs_net_charge_normalized": max(abs(row["net_charge_normalized"]) for row in charge_rows),
    }
    plot_summary = {
        "poisson_dx": [row["dx"] for row in dirichlet_rows],
        "dirichlet_phi_l2": [row["phi_l2_error"] for row in dirichlet_rows],
        "periodic_dx": [row["dx"] for row in periodic_rows],
        "periodic_phi_l2": [row["phi_l2_error"] for row in periodic_rows],
    }
    plot_path = write_plot(
        out_dir,
        dirichlet_plot,
        periodic_plot,
        reformulated_plot,
        gaussian_plot,
        profile_plot,
        plot_summary,
    )
    if plot_path is not None:
        summary["figure"] = str(plot_path)

    with (out_dir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    if args.assert_validation:
        assert_validation(summary)  # type: ignore[arg-type]

    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
