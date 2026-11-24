"""Calibrate and test the existing filtration solver against digitized column data."""

from dataclasses import asdict, replace
import csv
import json
from pathlib import Path

import numpy as np
from scipy.linalg.lapack import dgttrf, dgttrs
from scipy.optimize import differential_evolution, minimize, minimize_scalar

from solver import Parameters, simulate


OUT = Path(__file__).resolve().parent / "results"
LENGTH = .40
RADIUS = .0225
POROSITY = .37
FLOW = 145e-6 / 60
VELOCITY = FLOW / (np.pi * RADIUS**2)
DISPERSIVITY = .008
DENSITY = 2600.0
END = 83
PV_TIME = POROSITY * LENGTH / VELOCITY
CELLS = 160
STEPS_PER_PV = 40


class LapackTransport:
    """LU factorization of the same tridiagonal flux matrix used by solver.py."""

    def __init__(self, parameters, cells, dx, dt):
        d = parameters.dispersion / dx
        w = dt / dx
        v = parameters.velocity
        diagonal = np.full(cells, parameters.mobile_porosity + w*(v + 2*d))
        diagonal[0] += w*d
        diagonal[-1] -= w*d
        lower = np.full(cells-1, -w*(v+d))
        upper = np.full(cells-1, -w*d)
        self.inlet = w*(v+2*d)*parameters.inlet
        *self.factors, info = dgttrf(lower, diagonal, upper)
        if info:
            raise ValueError(f"LAPACK factorization failed: {info}")

    def solve(self, rhs):
        rhs[0] += self.inlet
        x, info = dgttrs(*self.factors, rhs[:, None], overwrite_b=1)
        if info:
            raise ValueError(f"LAPACK solve failed: {info}")
        return x[:, 0]


def base_parameters(feed, immobile, exchange):
    return Parameters(mobile_porosity=POROSITY-immobile, immobile_porosity=immobile,
                      velocity=VELOCITY, dispersion=DISPERSIVITY*VELOCITY,
                      inlet=feed/DENSITY, exchange=exchange)


def run_trace(parameters, query, *, cells=CELLS, steps_per_pv=STEPS_PER_PV, full=False):
    query = np.asarray(query)
    indexes = query * steps_per_pv
    steps = np.concatenate((np.floor(indexes), np.ceil(indexes))).astype(int)
    if full:
        steps = np.concatenate((steps, np.arange(2, END*steps_per_pv+1, 2),
                                [END*steps_per_pv]))
    steps = np.unique(steps[steps>0])
    dt = PV_TIME/steps_per_pv
    # Reuse the canonical transport and kinetic update, replacing only its linear solve.
    grid, states, budget = simulate(parameters, LENGTH, LENGTH/cells, dt,
                                    tuple(steps*dt), transport_type=LapackTransport, track_budget=full)
    time = np.r_[0, steps / steps_per_pv]
    ratio = np.r_[0, [states[t][0,-1]/parameters.inlet for t in states]]
    answer = np.interp(query, time, ratio)
    if full:
        terminal = END*steps_per_pv*dt
        return answer, time, ratio, grid, states[terminal], budget[terminal]
    return answer


def fit_parameters(rows, immobile, exchange):
    p = base_parameters(.25, immobile, exchange)
    early = [r for r in rows if r["source_figure"]=="2"]
    main = [r for r in rows if r["split"]=="fit"]

    def coords(data):
        return np.array([r["pore_volumes"] for r in data]), np.array([r["observed"] for r in data])

    x,y = coords(early)
    def loss_initial(beta):
        candidate = replace(p, deposition_initial=beta, threshold=1., capacity=2.)
        return np.mean((run_trace(candidate,x)-y)**2)

    initial = minimize_scalar(loss_initial, bounds=(.0005,.022), method="bounded",
                              options={"xatol":1e-7})
    p = replace(p, deposition_initial=float(initial.x))
    x,y = coords(main)
    threshold_min = 1.05*p.deposition_initial*p.inlet*3*PV_TIME
    bounds=np.log10(np.array([
        (.00035,.028), (.00008,.014), (threshold_min,.0045), (.00025,.023),
    ])).tolist()

    def decode(log_values):
        beta_a,beta_d,threshold,capacity_above_threshold = 10**np.asarray(log_values)
        return replace(p, deposition_late=float(beta_a), release=float(beta_d),
                       threshold=float(threshold), capacity=float(threshold+capacity_above_threshold))

    def loss(log_values):
        return np.mean((run_trace(decode(log_values),x)-y)**2)

    candidates = []
    search_rmse = []
    for seed in (20261124, 20261125):
        global_fit = differential_evolution(loss, bounds, maxiter=14, popsize=5, seed=seed,
                                             polish=False, tol=.001, atol=0)
        polished = minimize(loss, global_fit.x, method="Powell", bounds=bounds,
                            options={"maxfev":60, "xtol":.0003})
        candidates.extend((global_fit, polished))
        search_rmse.append(float(np.sqrt(min(global_fit.fun, polished.fun))))
    best = min(candidates, key=lambda result: result.fun)
    return decode(best.x), {
        "early_evaluations": int(initial.nfev),
        "later_evaluations": sum(int(result.nfev) for result in candidates),
        "search_rmse": search_rmse,
        "early_rmse": float(np.sqrt(initial.fun)),
        "calibration_rmse": float(np.sqrt(best.fun)),
        "optimized_log_bounds": bounds,
        "optimized_log_parameters": best.x.tolist(),
    }


def summarize(residual):
    errors=np.asarray(residual)
    return {
        "count":len(errors), "rmse":float(np.sqrt(np.mean(errors**2))),
        "bias":float(errors.mean()), "mae":float(np.mean(np.abs(errors))),
    }


def evaluate(rows, parameters, *, cells=CELLS, steps_per_pv=STEPS_PER_PV):
    by_feed={}
    all_predictions=[]
    for feed in (.25,.50):
        selection=[r for r in rows if r["feed_g_L"]==feed]
        btc=[r for r in selection if r["observable"]=="breakthrough"]
        dep=[r for r in selection if r["observable"]=="deposit"]
        p=replace(parameters,inlet=feed/DENSITY)
        bt,t,curve,x,state,budget=run_trace(p,[r["pore_volumes"] for r in btc],
                                            cells=cells,steps_per_pv=steps_per_pv,full=True)
        depth=np.array([r["distance_m"] for r in dep])
        rho=np.interp(depth,x,state[2])
        predictions=[{**row,"prediction":float(value)} for row,value in zip(btc+dep,np.r_[bt,rho])]
        all_predictions.extend(predictions)
        observed=np.array([r["observed"] for r in dep])
        metrics={
            "breakthrough_holdout":summarize([r["prediction"]-r["observed"] for r in predictions
                                             if r["observable"]=="breakthrough" and r["split"]=="holdout"]),
            "deposit":summarize(rho-observed),
            "observed_deposit_span_integral_m":float(np.trapezoid(observed,depth)),
            "predicted_deposit_span_integral_m":float(np.trapezoid(rho,depth)),
            "observed_deposit_first_last_ratio":float(observed[0]/observed[-1]),
            "predicted_deposit_first_last_ratio":float(rho[0]/rho[-1]),
            "mass_residual":budget["relative_mass_residual"],
        }
        if feed==.25:
            metrics["calibration"]=summarize([r["prediction"]-r["observed"] for r in predictions
                                              if r["split"]=="fit"])
        by_feed[f"{feed:.2f}"]={"metrics":metrics,"pv":t.tolist(),"breakthrough":curve.tolist(),
                                "x":x.tolist(),"deposit":state[2].tolist()}
    return by_feed,all_predictions


def check_implementation():
    p=replace(base_parameters(.25,.05,.005), threshold=.000012, capacity=.000065)
    times=(10.,250.,800.)
    defaults=simulate(p,LENGTH,.005,1.0,times)
    fast=simulate(p,LENGTH,.005,1.0,times,transport_type=LapackTransport)
    for t in defaults[1]:
        np.testing.assert_allclose(defaults[1][t],fast[1][t],rtol=5e-12,atol=1e-15)
    assert defaults[1][times[-1]][2].max() >= p.capacity-1e-13


def main():
    check_implementation()
    with (OUT/"column_observations.csv").open(encoding="utf-8",newline="") as f:
        rows=list(csv.DictReader(f))
    for r in rows:
        for field in ("feed_g_L","pore_volumes","distance_m","observed"):
            r[field]=float(r[field]) if r[field] else None
    assert len(rows)==89 and sum(r["split"]=="fit" for r in rows)==32
    report={}
    predictions=[]
    for name, immobile, exchange in (("exchange",.05,.005),("no_exchange",0.,0.)):
        print(f"Calibrating {name}",flush=True)
        p,fit=fit_parameters(rows,immobile,exchange)
        by_feed,obs=evaluate(rows,p)
        fine,_=evaluate(rows,p,cells=2*CELLS,steps_per_pv=2*STEPS_PER_PV)
        predictions.extend({"model":name,**r} for r in obs)
        report[name]={"parameters":asdict(p),"fit":fit,"cases":by_feed,
                      "fine_metrics":{k:v["metrics"] for k,v in fine.items()}}
        print(name,fit, {k:v["metrics"] for k,v in by_feed.items()},flush=True)
    metadata={"length_m":LENGTH,"radius_m":RADIUS,"porosity":POROSITY,"flow_m3_s":FLOW,
              "velocity_m_s":VELOCITY,"dispersivity_m":DISPERSIVITY,"density_kg_m3":DENSITY,
              "pore_volume_time_s":PV_TIME,"cells":CELLS,"steps_per_pv":STEPS_PER_PV,
              "end_pore_volumes":END,"models":report}
    (OUT/"column_validation.json").write_text(json.dumps(metadata,indent=2)+"\n",encoding="utf-8")
    with (OUT/"column_predictions.csv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(predictions[0]))
        writer.writeheader()
        writer.writerows(predictions)
    print("Saved column validation results.",flush=True)


if __name__=="__main__":
    main()
