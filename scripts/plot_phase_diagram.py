#!/usr/bin/env python
"""Publication-quality phase-diagram figure from a scan CSV.

Draws the stability classification as a filled map on the regular scan grid,
contours of the optimized ansatz action, the metastable/unstable boundary,
the benchmark point, and the PDG 2022 1-sigma/2-sigma experimental region.

Usage:
    python scripts/plot_phase_diagram.py results/numerical_data.csv \
        --output figures/phase_diagram.png [--region sm|full]
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse

# PDG 2022 central values and 1-sigma uncertainties (GeV).
PDG_MH, PDG_MH_ERR = 125.20, 0.11
PDG_MT, PDG_MT_ERR = 172.57, 0.29
# Benchmark point used throughout the code and reference tables.
BENCH_MH, BENCH_MT = 125.1, 173.1

SM_REGION = dict(mh=(110.0, 140.0), mt=(155.0, 185.0))

STATUS_COLORS = {1: "#2166ac", 2: "#f7f7f7", 3: "#b2182b", 4: "#878787"}
STATUS_LABELS = {1: "stable", 2: "metastable", 3: "unstable", 4: "non-perturbative"}


def regular_grid(df):
    """Reshape the Mt-major scan CSV into (Mt, Mh) value grids."""
    mts = np.sort(df["Mt"].unique())
    mhs = np.sort(df["Mh_calc"].unique())
    step_mh = np.min(np.diff(mhs)) if len(mhs) > 1 else 1.0
    shape = (len(mts), len(mhs))
    order = np.lexsort((df["Mh_calc"], df["Mt"]))
    grids = {}
    for col in ("Stability", "S_exact", "S_approx"):
        vals = df[col].to_numpy()[order]
        if vals.size != shape[0] * shape[1]:
            raise SystemExit(f"CSV is not a complete regular grid ({vals.size} vs {shape[0]*shape[1]})")
        grids[col] = vals.reshape(shape)
    return mts, mhs, grids


def plot(input_csv, output, region):
    df = pd.read_csv(input_csv)
    if region == "sm":
        sel = ((df["Mt"] >= SM_REGION["mt"][0]) & (df["Mt"] <= SM_REGION["mt"][1]) &
               (df["Mh_calc"] >= SM_REGION["mh"][0]) & (df["Mh_calc"] <= SM_REGION["mh"][1]))
        df = df[sel]
    mts, mhs, g = regular_grid(df)

    fig, ax = plt.subplots(figsize=(7.0, 5.4))

    # Filled classification map (midpoint colors between grid cells).
    status = g["Stability"]
    from matplotlib.colors import ListedColormap, BoundaryNorm
    cmap = ListedColormap([STATUS_COLORS[i] for i in range(1, 5)])
    norm = BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)
    mh_edges = np.append(mhs - (mhs[1] - mhs[0]) / 2, mhs[-1] + (mhs[1] - mhs[0]) / 2)
    mt_edges = np.append(mts - (mts[1] - mts[0]) / 2, mts[-1] + (mts[1] - mts[0]) / 2)
    ax.pcolormesh(mh_edges, mt_edges, status, cmap=cmap, norm=norm,
                  shading="flat", rasterized=True)

    # Contours of the optimized ansatz action (where defined).
    S = np.where(g["S_exact"] > 0, g["S_exact"], np.nan)
    levels = [450, 600, 800, 1200, 2000, 4000]
    cs = ax.contour(mhs, mts, S, levels=levels, colors="k", linewidths=0.6, alpha=0.55)
    ax.clabel(cs, fmt="S=%d", fontsize=7, inline=True)

    # Metastable/unstable boundary (status 3 vs 2).
    unb = np.where(g["Stability"] == 3, 1.0, 0.0)
    if 0 < unb.sum() < unb.size:
        ax.contour(mhs, mts, unb, levels=[0.5], colors="#b2182b", linewidths=2.0)

    # Benchmark point and PDG 2022 uncertainty ellipses.
    ax.plot(BENCH_MH, BENCH_MT, marker="*", color="black", markersize=14,
            linestyle="none", zorder=6, markeredgecolor="white", markeredgewidth=0.5,
            label="benchmark point (125.1, 173.1)")
    for factor, style in ((1, "--"), (2, ":")):
        ax.add_patch(Ellipse((PDG_MH, PDG_MT), width=2 * factor * PDG_MH_ERR,
                             height=2 * factor * PDG_MT_ERR, fill=False,
                             edgecolor="black", linestyle=style, linewidth=1.0, zorder=5))
    ax.plot([], [], linestyle="--", color="black", label="PDG 2022 $1\\sigma$")
    ax.plot([], [], linestyle=":", color="black", label="PDG 2022 $2\\sigma$")

    ax.set_xlabel(r"$M_h$ [GeV]")
    ax.set_ylabel(r"$M_t$ [GeV]")
    ax.set_xlim(mhs[0], mhs[-1])
    ax.set_ylim(mts[0], mts[-1])

    from matplotlib.patches import Patch
    handles = [Patch(facecolor=STATUS_COLORS[s], label=STATUS_LABELS[s]) for s in (1, 2, 3, 4)]
    line_handles, line_labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles + line_handles[1:], loc="upper left", framealpha=0.9,
              fontsize=8)

    fig.tight_layout()
    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    fig.savefig(output, dpi=300)
    print(f"Saved {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_csv")
    parser.add_argument("--output", default="figures/phase_diagram.png")
    parser.add_argument("--region", choices=["sm", "full"], default="sm")
    args = parser.parse_args()
    plot(args.input_csv, args.output, args.region)


if __name__ == "__main__":
    main()
