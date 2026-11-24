"""Compare uniform and depth-dependent capture against the same column data."""

import csv
import json

import numpy as np

from column_validation import OUT
from plot_column_validation import add_legend, create_panels, observations
from plot_results import save


def main():
    comparisons = []
    for name, model, label, style, color in (
        ("column_validation", "exchange", "Uniform capture", "--", "#9D450C"),
        ("straining_validation", "exchange", "Fitted straining", "-", "#006A86"),
        ("straining_validation", "full_straining", "Fully depth-dependent capture", "-.", "#007455"),
    ):
        report = json.loads((OUT / f"{name}.json").read_text(encoding="utf-8"))
        comparisons.append((report["models"][model]["cases"], label, style, color))
    with (OUT / "column_observations.csv").open(encoding="utf-8", newline="") as file:
        rows = list(csv.DictReader(file))

    for kind, x_key, x_scale, y_scale, x_label, y_label, figure in (
        ("breakthrough", "pv", 1, 1, "Injected pore volumes", "Outlet concentration ratio",
         "figure-9-straining-breakthrough"),
        ("deposit", "x", 100, 1000, "Distance from inlet", "Scaled deposit fraction",
         "figure-10-straining-deposit"),
    ):
        fig, axes = create_panels()
        if kind == "breakthrough":
            fig.set_figheight(2.60)
            fig.subplots_adjust(bottom=.41)
        upper = max(
            max(max(case[kind]) for cases, *_ in comparisons for case in cases.values()),
            max(float(row["observed"]) for row in rows if row["observable"] == kind),
        )
        for panel, (feed, ax) in enumerate(zip((.25, .50), axes)):
            ax.set_title(("(a) Calibration feed", "(b) Transfer feed")[panel],
                         loc="left", fontweight="bold")
            if kind == "breakthrough" and panel == 0:
                ax.axvspan(60, 83, facecolor="#edf2f4", edgecolor="none", zorder=0)
            for cases, label, style, color in comparisons:
                case = cases[f"{feed:.2f}"]
                ax.plot(np.asarray(case[x_key])*x_scale, np.asarray(case[kind])*y_scale,
                        style, color=color, lw=1.35, label=label)
            for split, face, label in (("fit", "#20262D", "Fit observations"),
                                       ("holdout", "white", "Held-out observations")):
                x, y = observations(rows, feed, kind, split)
                if len(x):
                    ax.plot(x*x_scale, y*y_scale, linestyle="", marker="o", ms=3.5,
                            mfc=face, mec="#20262D", mew=.85, label=label)
            ax.set(xlim=(0, 85 if kind == "breakthrough" else 41),
                   ylim=(0, 1.08*upper*y_scale), xlabel=x_label)
        axes[0].set_ylabel(y_label)
        add_legend(fig, axes)
        save(fig, figure)


if __name__ == "__main__":
    main()
