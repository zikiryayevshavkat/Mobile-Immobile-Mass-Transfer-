"""Render the numerical figures from the archived solver output."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, PathPatch
from matplotlib.path import Path as PlotPath
import numpy as np

from solver import Parameters


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
STYLES = (
    ("#172C3C", "-"),
    ("#0073A3", "--"),
    ("#A34800", ":"),
    ("#007455", "-."),
)
plt.rcParams.update({
    "font.size": 8,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.labelcolor": "#172C3C",
    "text.color": "#172C3C",
    "mathtext.default": "regular",
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def save(fig, name):
    fig.savefig(RESULTS / f"{name}.png", dpi=600)
    fig.savefig(RESULTS / f"{name}.pdf")
    plt.close(fig)


def plot_schematic():
    fig = plt.figure(figsize=(5.36, 1.91))
    ax = fig.add_axes((0, 0, 1, 1), xlim=(0, 386), ylim=(0, 137))
    ax.axis("off")
    boundary = PlotPath(
        [(17, 18), (17, 128), (321, 128), (316, 103), (311, 93), (319, 74),
         (335, 47), (329, 40), (321, 18), (17, 18)],
        [PlotPath.MOVETO, PlotPath.LINETO, PlotPath.LINETO, *([PlotPath.CURVE4] * 6), PlotPath.CLOSEPOLY],
    )
    ax.add_patch(PathPatch(boundary, facecolor="#F1F7FA", edgecolor="#20262D", lw=.85))
    islands = (
        [(42, 111), (46, 123), (75, 112), (94, 111), (129, 108), (173, 120), (181, 108),
         (195, 84), (175, 87), (153, 88), (112, 90), (85, 87), (65, 80),
         (55, 77), (43, 99), (42, 111)],
        [(197, 78), (203, 101), (231, 110), (252, 114), (274, 119), (275, 94), (276, 66),
         (278, 42), (266, 62), (252, 65), (237, 70), (195, 60), (197, 78)],
        [(79, 59), (96, 74), (133, 77), (155, 64), (145, 58), (179, 66), (194, 60),
         (212, 50), (226, 33), (217, 33), (200, 36), (178, 37), (154, 32),
         (124, 28), (84, 29), (72, 34), (66, 36), (88, 50), (83, 51),
         (79, 54), (77, 56), (79, 59)],
    )
    for island in islands:
        path = PlotPath(island, [PlotPath.MOVETO] + [PlotPath.CURVE4] * (len(island) - 1))
        ax.add_patch(PathPatch(path, facecolor="#C3D3DA", edgecolor="#354E59", lw=.85))
    arrows = (
        [(45, 60), (60, 70), (81, 77), (108, 81)],
        [(172, 80), (189, 81), (201, 95), (210, 108), (217, 121), (239, 120), (256, 120)],
        [(185, 69), (216, 57), (244, 41), (285, 58)],
    )
    for points in arrows:
        path = PlotPath(points, [PlotPath.MOVETO] + [PlotPath.CURVE4] * (len(points) - 1))
        ax.add_patch(FancyArrowPatch(path=path, arrowstyle="-|>", mutation_scale=9, color="#12627D"))
    ax.annotate("", (116, 23), (70, 23), arrowprops={"arrowstyle": "-|>", "color": "#12627D"})
    ax.annotate("", (371, 18), (17, 18), arrowprops={"arrowstyle": "-|>", "color": "#20262D"})
    ax.text(29, 36, "Mobile\nliquid", va="center", fontsize=8)
    ax.text(130, 46, "Immobile liquid", ha="center", va="center", fontsize=8)
    ax.text(17, 6, "Inlet", ha="left", va="center", fontsize=8)
    ax.text(371, 6, "Distance", ha="right", va="center", fontsize=8)
    save(fig, "figure-1-schematic")


def plot_profiles(x, fields, labels, name):
    fig, axes = plt.subplots(1, 3, figsize=(5.36, 2.20))
    fig.subplots_adjust(left=.09, right=.957, bottom=.29, top=.84, wspace=.43)
    titles = ("(a) Mobile", "(b) Immobile", "(c) Deposit")
    limits = (1.05, .052, .074)
    scales = (Parameters().inlet, 1, 1)
    for panel, (ax, title, limit, scale) in enumerate(zip(axes, titles, limits, scales)):
        for state, label, (color, linestyle) in zip(fields, labels, STYLES):
            ax.plot(x, state[panel] / scale, label=label, color=color, ls=linestyle, lw=1.3)
        ax.set_title(title, loc="left", fontweight="bold", pad=6)
        ax.set(xlim=(0, .45), ylim=(0, limit), xticks=(0, .15, .30, .45), xlabel="Distance")
        ax.grid(axis="y", color="#e1e7ec", lw=.65)
        ax.set_axisbelow(True)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(labels), frameon=False,
               handlelength=2.3, columnspacing=1.2)
    save(fig, name)


def plot_benchmark(data):
    x = data["benchmark_x"]
    inlet = Parameters().inlet
    exact = data["benchmark_exact"] / inlet
    numerical = data["benchmark_numerical"] / inlet
    fig, axes = plt.subplots(1, 2, figsize=(5.36, 2.20))
    fig.subplots_adjust(left=.10, right=.97, bottom=.23, top=.82, wspace=.43)
    axes[0].plot(x, exact, color=STYLES[0][0], label="Ogata–Banks", lw=1.6)
    axes[0].plot(x, numerical, "--", color=STYLES[2][0], label="Numerical", lw=1.2)
    axes[0].set(title="(a) Normalized concentration")
    axes[0].legend(frameon=False, loc="upper right", fontsize=7)
    axes[1].plot(x, 1000 * (numerical - exact), color=STYLES[1][0], lw=1.4)
    axes[1].axhline(0, color="#8695A0", ls=":", lw=.8)
    axes[1].set(title="(b) Scaled signed error")
    for ax in axes:
        ax.set(xlabel="Distance", xlim=(0, .55))
        ax.grid(axis="y", color="#e1e7ec", lw=.65)
        ax.set_axisbelow(True)
    save(fig, "figure-2-benchmark")


def main():
    plot_schematic()
    with np.load(RESULTS / "profiles.npz") as data:
        plot_benchmark(data)
        rates = data["rates"]
        fields = data["fields"]
        x = data["x"]
        times = data["times"]
        standard = np.flatnonzero(np.isclose(rates, Parameters().exchange))[0]
        plot_profiles(x, fields[standard], [f"{time:.0f}" for time in times], "figure-3-time")
        plot_profiles(x, fields[:, -1], [f"{k:.3f}" for k in rates], "figure-4-exchange")


if __name__ == "__main__":
    main()
