"""Redraw the column comparison using digitized points and our model predictions."""

import csv
import json

import numpy as np

from plot_results import plt, save
from column_validation import OUT


def create_panels():
    fig, axes = plt.subplots(1, 2, figsize=(5.36, 2.45))
    fig.subplots_adjust(left=.115, right=.965, top=.84, bottom=.32, wspace=.32)
    for ax in axes:
        ax.grid(axis="y", color="#e1e7ec", lw=.65)
        ax.set_axisbelow(True)
    return fig, axes


def observations(rows, feed, kind, split):
    rows = [r for r in rows if float(r["feed_g_L"])==feed and r["observable"]==kind and r["split"]==split]
    axis = "pore_volumes" if kind=="breakthrough" else "distance_m"
    return np.array([float(r[axis]) for r in rows]), np.array([float(r["observed"]) for r in rows])


def add_curves(ax, models, feed, x_key, y_key, x_scale=1, y_scale=1):
    for name, label, style, color in (
        ("exchange", "With exchange", "-", "#006A86"),
        ("no_exchange", "Zero exchange", "--", "#9D450C"),
    ):
        data = models[name]["cases"][f"{feed:.2f}"]
        ax.plot(np.array(data[x_key])*x_scale, np.array(data[y_key])*y_scale, style,
                color=color, label=label, lw=1.25)


def add_legend(fig, axes):
    entries = {}
    for ax in axes:
        handles, labels = ax.get_legend_handles_labels()
        entries.update(zip(labels,handles))
    fig.legend(list(entries.values()), list(entries), loc="lower center", ncol=2, frameon=False)


def main():
    models = json.loads((OUT/"column_validation.json").read_text(encoding="utf-8"))["models"]
    with (OUT/"column_observations.csv").open(encoding="utf-8",newline="") as file:
        data=list(csv.DictReader(file))
    fig, axes = create_panels()
    max_concentration = max(max(case["breakthrough"]) for model in models.values()
                            for case in model["cases"].values())
    for i, (feed, ax) in enumerate(zip((.25,.50), axes)):
        ax.set_title(("(a) Calibration feed","(b) Transfer feed")[i],loc="left",fontweight="bold")
        if i==0:
            ax.axvspan(60,83,facecolor="#edf2f4",edgecolor="none",zorder=0)
        add_curves(ax,models,feed,"pv","breakthrough")
        fit=observations(data,feed,"breakthrough","fit")
        if len(fit[0]):
            ax.plot(*fit,linestyle="",marker="o",ms=3.1,color="#20262D",label="Fit observations")
        held=observations(data,feed,"breakthrough","holdout")
        ax.plot(*held,linestyle="",marker="o",ms=3.7,mfc="white",mew=.9,
                mec="#20262D",color="#20262D",label="Held-out observations")
        ax.set(xlim=(0,85),ylim=(0,max(.65,1.05*max_concentration)),xlabel="Injected pore volumes")
    axes[0].set_ylabel("Outlet concentration ratio")
    add_legend(fig,axes)
    save(fig,"figure-7-column-breakthrough")

    fig, axes = create_panels()
    max_deposit = max(
        max(max(case["deposit"]) for model in models.values() for case in model["cases"].values())*1000,
        max(float(r["observed"]) for r in data if r["observable"]=="deposit")*1000,
    )
    for i,(feed, ax) in enumerate(zip((.25,.50), axes)):
        ax.set_title(("(a) Calibration feed","(b) Transfer feed")[i],loc="left",fontweight="bold")
        add_curves(ax,models,feed,"x","deposit",100,1000)
        x,y=observations(data,feed,"deposit","holdout")
        ax.plot(x*100,y*1000,linestyle="",marker="o",ms=4,mfc="white",mew=.9,
                mec="#20262D",color="#20262D",label="Held-out observations")
        ax.set(xlim=(0,41),ylim=(0,1.08*max_deposit),xlabel="Distance from inlet")
    axes[0].set_ylabel("Scaled deposit fraction")
    add_legend(fig,axes)
    save(fig,"figure-8-column-deposit")


if __name__=="__main__":
    main()
