#!/usr/bin/env python
"""Figure 11: the narrow structures along the stability boundary at the
highest available sampling resolution (0.05 GeV).

Over the classification fills (green/yellow/red), two narrow features are
drawn explicitly:

- brown dots: ansatz-breakdown points (classified unstable with S <= 0),
  tracing the stability boundary;
- blue squares: points where the strict conformal estimate and the
  RG-improved evaluation classify differently.

The figure answers: how wide are these artifact/disagreement strips really?
At 0.05 GeV sampling both strips are one to two grid cells wide, i.e.
physically ≲ 0.1 GeV; they are features of the fixed-profile ansatz and of
the threshold criterion, not extended physical regions.
"""

import argparse
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_style import (apply_style, BENCH_MH, BENCH_MT, STATUS_LABELS)
from plot_phase_diagram import regular_grid

BOUNDARY_WINDOW = dict(mh=(110.0, 140.0), mt=(162.0, 182.0))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", default="figures/11_boundary_strip_0p05GeV.png")
    args = parser.parse_args()

    apply_style()
    num = pd.read_csv("data/numerical_boundary_0p05GeV.csv")
    ana = pd.read_csv("data/analytical_boundary_0p05GeV.csv")
    key = ["Mt", "Mh_calc"]
    m = num[key + ["Stability", "S_trial", "lambda_min", "method_flag"]].merge(
        ana[key + ["Stability"]], on=key, how="inner", suffixes=("_rg", "_conf"))
    mts, mhs, g = regular_grid(m, ["Stability_rg", "S_trial", "Stability_conf",
                                   "lambda_min", "method_flag"])

    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    from figure_style import STATUS_COLORS
    cmap = ListedColormap([STATUS_COLORS[i] for i in (1, 2, 3)])
    ax.pcolormesh(mhs, mts, np.where(g["Stability_rg"] == 4, np.nan,
                                     g["Stability_rg"]),
                  cmap=ListedColormap([STATUS_COLORS[i] for i in (1, 2, 3)]),
                  norm=BoundaryNorm([0.5, 1.5, 2.5, 3.5], 3),
                  shading="auto", rasterized=True)

    # Stability boundary as the lambda_min = 0 level set of the 0.05 GeV field.
    lam = np.where(g["Stability_rg"] != 4, g["lambda_min"], np.nan)
    if np.isfinite(lam).any() and (lam > 0).any() and (lam < 0).any():
        ax.contour(mhs, mts, lam, levels=[0.0], colors="black", linewidths=1.2)

    breakdown = g["method_flag"] == 3
    diff = g["Stability_rg"] != g["Stability_conf"]
    MM, HH = np.meshgrid(mts, mhs, indexing="ij")
    if diff.any():
        ax.plot(HH[diff], MM[diff], linestyle="none", marker="s", markersize=1.8,
                color="#2166ac", label="classification difference (conformal vs RG)")
    if breakdown.any():
        ax.plot(HH[breakdown], MM[breakdown], linestyle="none", marker="o",
                markersize=1.8, color="#7f2704", label="ansatz breakdown ($S \\leq 0$)")

    nbd = int(breakdown.sum())
    ndiff = int(diff.sum())
    ax.plot(BENCH_MH, BENCH_MT, marker="*", color="black", markersize=12,
            linestyle="none", markeredgecolor="white", markeredgewidth=0.5, zorder=6)

    handles = [Patch(facecolor=STATUS_COLORS[s], label=STATUS_LABELS[s])
               for s in (1, 2, 3)]
    handles += [
        Line2D([], [], color="black", lw=1.2,
               label=r"stability boundary ($\lambda_{\min}=0$)"),
        Line2D([], [], linestyle="none", marker="s", markersize=5, color="#2166ac",
               label="classification difference"),
        Line2D([], [], linestyle="none", marker="o", markersize=5, color="#7f2704",
               label="ansatz breakdown ($S \\leq 0$)"),
        Line2D([], [], linestyle="none", marker="*", markersize=10,
               markerfacecolor="black", markeredgecolor="white",
               label="benchmark (125.1, 173.1) GeV"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=7.5, framealpha=0.9,
              frameon=True, edgecolor="0.8")
    ax.text(0.98, 0.02,
            f"{nbd} breakdown, {ndiff} difference points "
            f"({ndiff / g['Stability_rg'].size:.2%} of the strip)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8)

    ax.set_xlabel(r"$M_h$ [GeV]")
    ax.set_ylabel(r"$M_t$ [GeV]")
    ax.set_xlim(mhs[0], mhs[-1])
    ax.set_ylim(mts[0], mts[-1])

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
