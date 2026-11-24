"""Factorial sensitivity cases for immobile storage and the deposition threshold."""

from dataclasses import replace
import csv
import json
from pathlib import Path

import numpy as np

from solver import Parameters, simulate


OUT = Path(__file__).resolve().parent / "results"
BASE = Parameters()
TIMES = (450, 900, 1350, 1800)
LENGTH = 0.60
DX, DT = 0.001, 0.25
POROSITIES = (0.10, 0.20, 0.30)
RATES = (0.0, 0.002, 0.005, 0.020)
THRESHOLDS = (0.003, 0.007, 0.015)
POINTS = (0.05, 0.10, 0.20)


def key(porosity, exchange, threshold):
    return f"{porosity:.2f}/{exchange:.3f}/{threshold:.3f}"


def summarize(grid, state, budget):
    c, n, rho = state
    storage = budget["storage"]
    deposit = rho.sum()
    return {
        "mobile_volume": storage[0],
        "immobile_volume": storage[1],
        "deposit_volume": storage[2],
        "immobile_share": storage[1] / sum(storage),
        "deposit_centroid_m": float(grid @ rho / deposit),
        "inlet_volume": budget["inlet"],
        "outlet_volume": budget["outlet"],
        "mass_residual": budget["relative_mass_residual"],
        "points": {
            f"{distance:.2f}": [float(np.interp(distance, grid, field)) for field in (c, n, rho)]
            for distance in POINTS
        },
    }


def run():
    OUT.mkdir(exist_ok=True)
    cases = {}
    profiles = {}

    def record(porosity, exchange, threshold):
        name = key(porosity, exchange, threshold)
        if name not in cases:
            p = replace(BASE, immobile_porosity=porosity, exchange=exchange, threshold=threshold)
            x, states, budgets = simulate(p, LENGTH, DX, DT, TIMES)
            cases[name] = {str(t): summarize(x, states[t], budgets[t]) for t in TIMES}
            profiles[name] = np.stack([states[t] for t in TIMES])
            print(f"Finished {name}", flush=True)
        return name

    for porosity in POROSITIES:
        for exchange in RATES:
            record(porosity, exchange, BASE.threshold)
    for threshold in THRESHOLDS:
        for exchange in (0, BASE.exchange):
            record(BASE.immobile_porosity, exchange, threshold)

    # The zero-exchange result is independent of immobile porosity.
    assert all(np.array_equal(profiles[key(p, 0, BASE.threshold)],
                              profiles[key(BASE.immobile_porosity, 0, BASE.threshold)])
               for p in POROSITIES)
    with np.load(OUT / "profiles.npz") as previous:
        grid = previous["x"]
        for exchange, archived in zip(previous["rates"], previous["fields"]):
            np.testing.assert_allclose(profiles[key(BASE.immobile_porosity, exchange, BASE.threshold)],
                                       archived, rtol=1e-11, atol=1e-13)
    refinements = {}
    for porosity, exchange, threshold in (
        (0.10, 0.005, 0.007), (0.30, 0.005, 0.007),
        (0.20, 0.005, 0.003), (0.20, 0.005, 0.015),
    ):
        p = replace(BASE, immobile_porosity=porosity, exchange=exchange, threshold=threshold)
        x, states, budgets = simulate(p, LENGTH, DX / 2, DT / 2, (TIMES[-1],))
        name = key(porosity, exchange, threshold)
        sample = summarize(x, states[TIMES[-1]], budgets[TIMES[-1]])
        coarse = cases[name][str(TIMES[-1])]
        refinements[name] = {
            "relative_deposit_change": abs(coarse["deposit_volume"] / sample["deposit_volume"] - 1),
            "relative_centroid_change": abs(coarse["deposit_centroid_m"] / sample["deposit_centroid_m"] - 1),
            "fine": sample,
        }
        print(f"Refined {name}", flush=True)

    assert max(cases[name]["1800"]["mass_residual"] for name in cases) < 1e-8
    assert max(row["relative_deposit_change"] for row in refinements.values()) < 0.01
    metadata = {
        "dx": DX, "dt": DT, "length": LENGTH, "times": TIMES,
        "case_key_order": ["immobile_porosity", "exchange_per_s", "deposit_threshold"],
        "cases": cases, "refinements": refinements,
    }
    (OUT / "sensitivity.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    np.savez_compressed(OUT / "sensitivity_profiles.npz", x=grid, times=TIMES,
                        cases=list(profiles), fields=np.stack(list(profiles.values())))
    with (OUT / "sensitivity_budgets.csv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(("immobile_porosity", "exchange_per_s", "deposit_threshold", "time_s",
                         "mobile_volume_per_area_m", "immobile_volume_per_area_m",
                         "deposit_volume_per_area_m", "immobile_share", "deposit_centroid_m",
                         "inlet_volume_per_area_m", "outlet_volume_per_area_m", "mass_residual"))
        for name, samples in cases.items():
            for t, metrics in samples.items():
                writer.writerow((*name.split("/"), t, *[v for k, v in metrics.items() if k != "points"]))
    print("Finished sensitivity cases and wrote source data.", flush=True)


if __name__ == "__main__":
    run()
