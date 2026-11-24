"""Conservative one-dimensional suspension filtration and verification cases."""

from dataclasses import dataclass, replace
import csv
import json
from math import erfc
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Parameters:
    mobile_porosity: float = 0.30
    immobile_porosity: float = 0.20
    velocity: float = 1.0e-4
    dispersion: float = 2.5e-6
    inlet: float = 0.05
    exchange: float = 0.005
    deposition_initial: float = 0.003
    deposition_late: float = 0.00777
    release: float = 0.005
    threshold: float = 0.007
    capacity: float = 0.070
    straining_fraction: float = 0.0
    grain_diameter: float = 0.000410
    straining_exponent: float = 0.432


def capture_weight(x, parameters):
    """Uniform and depth-dependent shares of the local capture coefficient."""
    fraction = parameters.straining_fraction
    if not 0 <= fraction <= 1:
        raise ValueError("The inlet straining fraction must be between zero and one.")
    if parameters.grain_diameter <= 0 or parameters.straining_exponent <= 0:
        raise ValueError("Straining grain diameter and exponent must be positive.")
    return (1 - fraction) + fraction*(1 + x/parameters.grain_diameter)**(-parameters.straining_exponent)


class TridiagonalTransport:
    """Factor the backward-Euler upwind advection and centered dispersion operator."""

    def __init__(self, parameters, cells, dx, dt):
        diffusion = parameters.dispersion / dx
        ratio = dt / dx
        v = parameters.velocity
        m = parameters.mobile_porosity
        self.lower = -ratio * (v + diffusion)
        self.upper = -ratio * diffusion
        self.inlet_term = ratio * (v + 2 * diffusion) * parameters.inlet
        diagonal = np.full(cells, m + ratio * (v + 2 * diffusion))
        diagonal[0] += ratio * diffusion
        diagonal[-1] -= ratio * diffusion
        self.multiplier = np.empty(cells)
        for i in range(1, cells):
            self.multiplier[i] = self.lower / diagonal[i - 1]
            diagonal[i] -= self.multiplier[i] * self.upper
        self.inverse_diagonal = 1 / diagonal

    def solve(self, rhs):
        rhs[0] += self.inlet_term
        for i in range(1, rhs.size):
            rhs[i] -= self.multiplier[i] * rhs[i - 1]
        rhs[-1] *= self.inverse_diagonal[-1]
        for i in range(rhs.size - 2, -1, -1):
            rhs[i] = (rhs[i] - self.upper * rhs[i + 1]) * self.inverse_diagonal[i]
        return rhs


def simulate(p, length, dx, dt, times, *, transport_type=TridiagonalTransport, track_budget=True):
    cells = round(length / dx)
    steps = {round(time / dt): time for time in times}
    if cells < 2 or not np.isclose(cells * dx, length):
        raise ValueError("Choose a length that contains an integer number of cells.")
    if not times or len(steps) != len(times) or min(steps) < 1 or any(
        not np.isclose(step * dt, time) for step, time in steps.items()
    ):
        raise ValueError("Snapshot times must be distinct positive multiples of dt.")
    grid = (np.arange(cells) + 0.5) * dx
    weight = capture_weight(grid, p)
    sink = p.immobile_porosity * p.exchange + weight.max()*max(p.deposition_initial, p.deposition_late)
    if dt * sink > p.mobile_porosity or dt * max(p.exchange, p.release) > 1:
        raise ValueError("The explicit kinetic update requires a smaller time step.")

    c = np.zeros(cells)
    n = np.zeros(cells)
    rho = np.zeros(cells)
    transport = transport_type(p, cells, dx, dt)
    snapshots, budgets = {}, {}
    incoming, outgoing, max_residual = 0.0, 0.0, 0.0
    for step in range(1, max(steps) + 1):
        delta_n = dt * p.exchange * (c - n)
        available = c * weight if p.straining_fraction else c
        deposition = np.where(
            rho <= p.threshold,
            p.deposition_initial * available,
            p.deposition_late * available - p.release * rho,
        )
        next_rho = np.clip(rho + dt * deposition, 0, p.capacity)
        # Debit the mobile phase by the actual storage changes after projection.
        rhs = p.mobile_porosity * c - p.immobile_porosity * delta_n - (next_rho - rho)
        c = transport.solve(rhs)
        n += delta_n
        rho = next_rho
        if track_budget:
            incoming += dt * (
                p.velocity * p.inlet + 2 * p.dispersion * (p.inlet - c[0]) / dx
            )
            outgoing += dt * p.velocity * c[-1]
            storage = dx * np.array([p.mobile_porosity * c.sum(), p.immobile_porosity * n.sum(), rho.sum()])
            residual = abs(storage.sum() - incoming + outgoing) / incoming
            max_residual = max(max_residual, residual)
        if step in steps:
            assert all(np.isfinite(field).all() and field.min() >= -1e-14 for field in (c, n, rho))
            snapshots[steps[step]] = np.stack((c, n, rho))
            if track_budget:
                budgets[steps[step]] = {
                    "storage": storage.tolist(), "inlet": incoming, "outlet": outgoing,
                    "relative_mass_residual": max_residual,
                }
    return grid, snapshots, budgets


def ogata_banks(x, t, parameters):
    """Semi-infinite, initially clean, constant-concentration inlet solution."""
    diffusivity = parameters.dispersion / parameters.mobile_porosity
    speed = parameters.velocity / parameters.mobile_porosity
    scale = 2 * np.sqrt(diffusivity * t)
    erfc_values = np.vectorize(erfc, otypes=[float])
    return 0.5 * parameters.inlet * (
        erfc_values((x - speed * t) / scale)
        + np.exp(speed * x / diffusivity) * erfc_values((x + speed * t) / scale)
    )


def nested_error(coarse, fine, time):
    """Relative mean absolute errors on cell averages of the finer solution."""
    _, a, _ = coarse
    _, b, _ = fine
    shape = a[time].shape
    ratio = b[time].shape[1] // shape[1]
    target = b[time].reshape(*shape, ratio).mean(axis=2)
    return (np.abs(a[time] - target).sum(axis=1) / np.abs(target).sum(axis=1)).tolist()


def run(out):
    out.mkdir(parents=True, exist_ok=True)
    base = Parameters()
    times = (450, 900, 1350, 1800)
    rates = (0.0, 0.002, 0.005, 0.020)
    length, dx, dt = 0.60, 0.001, 0.25
    end = times[-1]
    profiles, sweep = [], []
    for k in rates:
        grid, states, budget = simulate(replace(base, exchange=k), length, dx, dt, times)
        profiles.append(np.stack([states[t] for t in times]))
        sample = [float(np.interp(0.10, grid, field)) for field in states[end]]
        sweep.append({"k": k, "at_0_10_m": sample, **budget[end]})
        if k == base.exchange:
            production = grid, states, budget
        print(f"Finished coefficient k={k:.3f}", flush=True)

    # Isolate one discretization variable at a time in the coupled problem.
    coarse_space = simulate(base, length, 0.004, dt, (end,))
    medium_space = simulate(base, length, 0.002, dt, (end,))
    coarse_time = simulate(base, length, dx, 1.00, (end,))
    medium_time = simulate(base, length, dx, 0.50, (end,))
    spatial = [nested_error(coarse_space, medium_space, end), nested_error(medium_space, production, end)]
    temporal = [nested_error(coarse_time, medium_time, end), nested_error(medium_time, production, end)]
    far_grid, far_states, _ = simulate(replace(base, exchange=0), 0.90, dx, dt, (end,))
    # The zero-exchange case is used here because its mobile front travels farthest.
    window = grid <= 0.45
    field_scale = np.array([base.inlet, base.inlet, base.capacity])
    domain_error = (np.max(np.abs(
        profiles[0][-1][:, window] - far_states[end][:, :len(grid)][:, window]
    ), axis=1) / field_scale).tolist()
    assert np.allclose(far_grid[:len(grid)], grid, rtol=0, atol=1e-14)

    bench = replace(base, exchange=0, deposition_initial=0, deposition_late=0, release=0)
    benchmark_rows = []
    for h, tau in ((0.008, 1.0), (0.004, 0.5), (0.002, 0.25), (0.001, 0.125)):
        _, values, _ = simulate(bench, 1.20, h, tau, (600,))
        x = (np.arange(values[600].shape[1]) + 0.5) * h
        # Three-point quadrature averages the analytical solution over each cell.
        offsets = np.sqrt(3 / 5) * h / 2
        values_exact = (
            5 * ogata_banks(x - offsets, 600, bench)
            + 8 * ogata_banks(x, 600, bench)
            + 5 * ogata_banks(x + offsets, 600, bench)
        ) / 18
        numerical = values[600][0]
        rel_l1 = np.abs(numerical - values_exact).sum() / np.abs(values_exact).sum()
        benchmark_rows.append({"dx": h, "dt": tau, "relative_l1": float(rel_l1)})
    orders = np.log2(np.asarray([row["relative_l1"] for row in benchmark_rows[:-1]]) /
                     np.asarray([row["relative_l1"] for row in benchmark_rows[1:]]))
    assert min(orders) > 0.8
    assert max(row["relative_mass_residual"] for row in sweep) < 1e-8
    assert np.all(profiles[0][:, 1, :] == 0)

    np.savez_compressed(out / "profiles.npz", x=grid, times=times, rates=rates, fields=profiles,
                        benchmark_x=x, benchmark_exact=values_exact, benchmark_numerical=numerical)
    (out / "verification.json").write_text(json.dumps({
        "sweep": sweep,
        "spatial_relative_l1": spatial,
        "temporal_relative_l1": temporal,
        "domain_normalized_linf": domain_error,
        "benchmark": benchmark_rows,
        "benchmark_orders": orders.tolist(),
    }, indent=2) + "\n", encoding="utf-8")
    with (out / "profiles.csv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(("k_per_s", "time_s", "x_m", "mobile_c", "immobile_n", "deposit_rho"))
        for k, case in zip(rates, profiles):
            for time, state in zip(times, case):
                writer.writerows((k, time, x_pos, *values) for x_pos, values in zip(grid, state.T))
    print("Finished verification and wrote profile data.", flush=True)


if __name__ == "__main__":
    run(Path(__file__).resolve().parent / "results")
