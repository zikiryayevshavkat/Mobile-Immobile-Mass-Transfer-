"""Calibrate a Bradford-type capture modifier on the original outlet-only split."""

from argparse import ArgumentParser
import csv
from dataclasses import asdict, replace
import json

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import differential_evolution, minimize, minimize_scalar

from column_validation import (
    CELLS, LENGTH, OUT, PV_TIME, STEPS_PER_PV, base_parameters, evaluate, run_trace,
)
from solver import Parameters, capture_weight


SEEDS = (20261126, 20261127)
MODELS = {"exchange": (.05, .005), "no_exchange": (0., 0.), "full_straining": (.05, .005)}
BASE_REPORT = OUT / "column_validation.json"
FINAL_REPORT = OUT / "straining_validation.json"


def observations():
    with (OUT / "column_observations.csv").open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        for field in ("feed_g_L", "pore_volumes", "distance_m", "observed"):
            row[field] = float(row[field]) if row[field] else None
    assert len(rows) == 89 and sum(r["split"] == "fit" for r in rows) == 32
    return rows


def coords(data):
    return np.array([r["pore_volumes"] for r in data]), np.array([r["observed"] for r in data])


def fit(rows, name, original, *, fixed_fraction=None):
    immobile, exchange = MODELS[name]
    p = base_parameters(.25, immobile, exchange)
    grid = (np.arange(CELLS)+.5)*LENGTH/CELLS
    t_early, y_early = coords([r for r in rows if r["source_figure"] == "2"])
    t_fit, y_fit = coords([r for r in rows if r["split"] == "fit"])
    uniform = Parameters(**original["parameters"])
    shape = capture_weight(grid, replace(p, straining_fraction=1))
    shape_mean = float(shape.mean())
    shape_max = float(shape.max())

    def weights(fraction):
        return (1-fraction+fraction*shape_mean, 1-fraction+fraction*shape_max)

    def early_rate(fraction):
        average_weight, _ = weights(fraction)
        candidate = replace(p, straining_fraction=fraction, threshold=1., capacity=2.)
        def error(average_beta):
            stage = replace(candidate, deposition_initial=average_beta/average_weight)
            return np.mean((run_trace(stage, t_early) - y_early)**2)
        result = minimize_scalar(error, bounds=(.0005, .022), method="bounded",
                                 options={"xatol":1e-8})
        return float(result.x)

    fractions = np.array([0,.20,.40,.60,.75,.85,.90,.94,.96,.98,.99,1])
    if fixed_fraction is None:
        early = PchipInterpolator(fractions, [early_rate(x) for x in fractions])
    else:
        average_initial_at_fixed = early_rate(fixed_fraction)
        early = lambda _: average_initial_at_fixed

    bounds = [(-3.4559319556497243,-1.5528419686577808),
              (-4.096910013008056,-1.853871964321762),
              (0,1), (-3.6020599913279625,-1.6382721639824072)]
    if fixed_fraction is None:
        bounds.append((0,1))

    def decode(x):
        fraction = float(x[-1]) if fixed_fraction is None else fixed_fraction
        average_weight, largest_weight = weights(fraction)
        beta_r = float(early(fraction))/average_weight
        threshold_min = 1.05*beta_r*largest_weight*p.inlet*3*PV_TIME
        threshold = threshold_min*(.0045/threshold_min)**float(x[2])
        return replace(p, deposition_initial=beta_r, deposition_late=10**float(x[0])/average_weight,
                       release=10**float(x[1]), threshold=threshold,
                       capacity=threshold+10**float(x[3]), straining_fraction=fraction)

    def loss(x):
        candidate = decode(x)
        # Out-of-range sink values can occur at full straining with the largest rate bound.
        dt = PV_TIME/STEPS_PER_PV
        max_weight = weights(candidate.straining_fraction)[1]
        sink = max_weight*max(candidate.deposition_initial, candidate.deposition_late)
        if dt*(candidate.immobile_porosity*candidate.exchange+sink)>candidate.mobile_porosity:
            return 1.
        return float(np.mean((run_trace(candidate,t_fit)-y_fit)**2))

    candidates = []
    if fixed_fraction is None:
        minimum = 1.05*float(early(0))*p.inlet*3*PV_TIME
        prior = [np.log10(uniform.deposition_late), np.log10(uniform.release),
                 np.clip(np.log(uniform.threshold/minimum)/np.log(.0045/minimum),0,1),
                 np.log10(uniform.capacity-uniform.threshold), 0.]
        candidates.append((loss(prior),prior))
    search = []
    for seed in SEEDS:
        fit_result = differential_evolution(loss,bounds,maxiter=15,popsize=5,seed=seed,polish=False,
                                            atol=0,tol=.001,x0=candidates[0][1] if candidates else None)
        local = minimize(loss,fit_result.x,method="Powell",bounds=bounds,
                         options={"maxfev":70,"xtol":.0005})
        best = min((fit_result,local),key=lambda trial:trial.fun)
        candidates.append((float(best.fun),best.x.tolist()))
        search.append(float(np.sqrt(best.fun)))
    mse, best = min(candidates,key=lambda item:item[0])
    return decode(best), {"rmse":float(np.sqrt(mse)), "search_rmse":search,"bounds":bounds,"solution":best}


def run_model(name):
    rows = observations()
    original_name = "exchange" if name == "full_straining" else name
    original = json.loads(BASE_REPORT.read_text(encoding="utf-8"))["models"][original_name]
    print(f"Fitting straining variant with {name} pore storage",flush=True)
    fixed = 1. if name == "full_straining" else None
    p,fit_info = fit(rows,name,original,fixed_fraction=fixed)
    by_feed,predictions = evaluate(rows,p)
    fine,_ = evaluate(rows,p,cells=2*CELLS,steps_per_pv=2*STEPS_PER_PV)
    result = {"parameters":asdict(p),"fit":fit_info,"cases":by_feed,
              "fine_metrics":{k:v["metrics"] for k,v in fine.items()}, "predictions":predictions}
    print(name, asdict(p), fit_info, {k:v["metrics"] for k,v in by_feed.items()},flush=True)
    return result


def write_report(models):
    old = json.loads(BASE_REPORT.read_text(encoding="utf-8"))
    old["models"] = models
    old["early_fit_method"] = "PCHIP interpolation of initial-rate fits across straining fraction"
    old["straining_shape"] = "Bradford-type depth factor, fixed exponent 0.432 and grain diameter 0.000410 m"
    old["full_straining_case"] = "Prescribed exchange with the inlet straining fraction fixed to one as a structural diagnostic"
    rows = []
    for name, result in models.items():
        rows.extend({"model":name,**row} for row in result.pop("predictions"))
    FINAL_REPORT.write_text(json.dumps(old,indent=2)+"\n",encoding="utf-8")
    with (OUT/"straining_predictions.csv").open("w",encoding="utf-8",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--model",choices=MODELS)
    parser.add_argument("--combine",action="store_true")
    args=parser.parse_args()
    if args.model:
        result=run_model(args.model)
        (OUT/f".straining-{args.model}.json").write_text(json.dumps(result),encoding="utf-8")
    elif args.combine:
        models={}
        for name in MODELS:
            partial=OUT/f".straining-{name}.json"
            models[name]=json.loads(partial.read_text(encoding="utf-8"))
        write_report(models)
        for name in MODELS:
            (OUT/f".straining-{name}.json").unlink()
    else:
        write_report({name: run_model(name) for name in MODELS})


if __name__ == "__main__":
    main()
