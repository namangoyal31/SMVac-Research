#!/usr/bin/env python
"""Figure 07: where the calculation is reliable.

Validity map over the full (Mh, Mt) plane, distinguishing (see
docs/limitations.md):

    white         : absolutely stable - no action is computed (none needed)
    light gray    : genuine ansatz action available (0 < S < 4.9e5)
    light orange  : kinetic-shortcut region (|lambda_R| <~ 3.3e-4): the
                    action is returned without the potential integral and is
                    only a rough overestimate
    dark orange   : fence points (no radius in the search bracket gives
                    lambda_R < 0; the exported S is the 1e100 sentinel)
    brown dots    : ansatz breakdown (classified unstable with S <= 0)
    black         : non-perturbative (|lambda| or y_t > 4 pi before the
                    Planck scale): the RGEs themselves leave their range of
                    validity

All categories are computed from the committed full-plane dataset; no point
is discarded silently.
"""

import argparse
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_style import apply_style, BENCH_MH, BENCH_MT
from plot_phase_diagram import regular_grid

# Category codes on the grid.
CAT_NO_ACTION, CAT_RELIABLE, CAT_SHORTCUT, CAT_FENCE, CAT_NONPERT = 0, 1, 2, 3, 4


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", default="figures/07_ansatz_breakdown.png")
    args = parser.parse_args()

    apply_style()
    df = pd.read_csv("data/numerical_full_plane.csv")
    mts, mhs, g = regular_grid(df, ["Stability", "S_exact", "S_kinetic"])

    status = g["Stability"]
    cat = np.full_like(status, CAT_RELIABLE, dtype=float)
    cat[status == 1] = CAT_NO_ACTION
    shortcut = (status != 1) & (g["S_exact"] > 0) & (g["S_exact"] >= 4.9e5) & (g["S_kinetic"] < 1e6)
    fence = (status != 1) & ((g["S_exact"] >= 1e99) | (g["S_kinetic"] >= 1e6))
    cat[shortcut] = CAT_SHORTCUT
    cat[fence] = CAT_FENCE
    cat[status == 4] = CAT_NONPERT
    breakdown = (status == 3) & (g["S_exact"] <= 0)
    cat[breakdown] = CAT_NO_ACTION  # fill under the dots (they sit in the unstable band)

    fig, ax = plt.subplots(figsize=(7.4, 6.4))
    cmap = ListedColormap(["white", "0.78", "#fdae6b", "#d94801", "black"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)
    ax.pcolormesh(mhs, mts, cat, cmap=cmap, norm=norm, shading="auto",
                  rasterized=True)

    MM, HH = np.meshgrid(mts, mhs, indexing="ij")
    if breakdown.any():
        ax.plot(HH[breakdown], MM[breakdown], linestyle="none", marker="o",
                markersize=2.4, color="#7f2704", label="ansatz breakdown ($S \\leq 0$)")

    ax.plot(BENCH_MH, BENCH_MT, marker="*", color="black", markersize=13,
            linestyle="none", markeredgecolor="white", markeredgewidth=0.5,
            zorder=6)

    handles = [
        Patch(facecolor="white", edgecolor="0.8", label="absolutely stable (no action computed)"),
        Patch(facecolor="0.78", label="genuine ansatz action available"),
        Patch(facecolor="#fdae6b", label="kinetic shortcut ($|\\lambda_R| \\lesssim 3\\times10^{-4}$)"),
        Patch(facecolor="#d94801", label="fence: no valid action in search bracket"),
        Line2D([], [], linestyle="none", marker="o", markersize=5, color="#7f2704",
               label="ansatz breakdown ($S \\leq 0$)"),
        Patch(facecolor="black", label="non-perturbative (RGEs out of range)"),
        Line2D([], [], linestyle="none", marker="*", markersize=10,
               markerfacecolor="black", markeredgecolor="white",
               label="benchmark (125.1, 173.1) GeV"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=7.0, framealpha=0.95,
              frameon=True, edgecolor="0.8")

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
