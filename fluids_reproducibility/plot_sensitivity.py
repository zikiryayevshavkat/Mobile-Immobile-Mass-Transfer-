"""Plot the two sensitivity figures using the solver's saved budgets."""

import json

from plot_results import STYLES, plt, save
from sensitivity import BASE, OUT, POROSITIES, RATES, THRESHOLDS, key


def panels():
    figure, axes = plt.subplots(1, 2, figsize=(5.36, 2.28))
    figure.subplots_adjust(left=.11, right=.961, bottom=.31, top=.83, wspace=.35)
    for ax in axes:
        ax.grid(axis="y", color="#e1e7ec", lw=.65)
        ax.set_axisbelow(True)
    return figure, axes


def finish(figure, axes, name, xlabel, ticks, fmt):
    for ax in axes:
        ax.set(xlabel=xlabel, xticks=ticks, xticklabels=[fmt.format(v) for v in ticks])
        ax.margins(x=.12, y=.15)
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=len(labels), frameon=False,
                  handlelength=2.2, columnspacing=1.4)
    save(figure, name)


def main():
    cases = json.loads((OUT / "sensitivity.json").read_text(encoding="utf-8"))["cases"]

    def sample(porosity, exchange, threshold):
        return cases[key(porosity, exchange, threshold)]["1800"]

    fig, axes = panels()
    baseline = sample(BASE.immobile_porosity, 0, BASE.threshold)["deposit_volume"]
    for rate, (color, line), marker in zip(RATES[1:], STYLES[1:], ("o", "s", "^")):
        values = [sample(p, rate, BASE.threshold) for p in POROSITIES]
        axes[0].plot(POROSITIES, [100 * (1 - a["deposit_volume"] / baseline) for a in values],
                     marker=marker, ms=4, color=color, ls=line, label=f"{rate:.3f}")
        axes[1].plot(POROSITIES, [100 * a["immobile_share"] for a in values],
                     marker=marker, ms=4, color=color, ls=line)
    axes[0].set(title="(a) Deposit reduction", ylabel="Percent")
    axes[1].set(title="(b) Immobile storage share", ylabel="Percent")
    finish(fig, axes, "figure-5-storage-sensitivity", "Immobile pore fraction", POROSITIES, "{:.2f}")

    fig, axes = panels()
    for rate, (color, line), marker in zip((0, BASE.exchange), STYLES, ("o", "s")):
        values = [sample(BASE.immobile_porosity, rate, threshold) for threshold in THRESHOLDS]
        axes[0].plot(THRESHOLDS, [a["deposit_volume"] * 1000 for a in values],
                     marker=marker, ms=4, color=color, ls=line, label=f"{rate:.3f}")
        axes[1].plot(THRESHOLDS, [a["deposit_centroid_m"] * 100 for a in values],
                     marker=marker, ms=4, color=color, ls=line)
    axes[0].set(title="(a) Integrated deposit", ylabel="Equivalent thickness")
    axes[1].set(title="(b) Deposit centroid", ylabel="Distance from inlet")
    finish(fig, axes, "figure-6-threshold-sensitivity", "Deposition threshold", THRESHOLDS, "{:.3f}")


if __name__ == "__main__":
    main()
