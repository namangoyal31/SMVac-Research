#!/usr/bin/env python
"""Figure 05: classification difference map — strict conformal estimate vs
RG-improved ansatz.

For every grid point of the full (Mh, Mt) plane, the vacuum classification
from the strict conformal estimate S_approx = 8 pi^2/(3|lambda_min|) is
compared with the classification from the RG-improved ansatz action
S_exact. Only points where the two classifications differ are highlighted:

    blue  : RG-improved classifies MORE stable than the conformal estimate
    red   : RG-improved classifies LESS stable
    black : either point is an ansatz breakdown / out-of-validity point

This is a comparison map, not a statement that either method is exact: the
RG-improved evaluation incorporates the running of the couplings across the
bubble profile, which the strict conformal limit neglects. Differences
concentrate along the classification boundaries, where the two actions
straddle the metastability threshold.
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
from figure_style import apply_style, BENCH_MH, BENCH_MT
from plot_phase_diagram import regular_grid


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", default="figures/05_conformal_vs_rg_difference.png")
    args = parser.parse_args()

    apply_style()
    num = pd.read_csv("data/numerical_full_plane.csv")
    ana = pd.read_csv("data/analytical_full_plane.csv")

    key = ["Mt", "Mh_calc"]
    merged = num[key + ["Stability", "S_exact"]].merge(
        ana[key + ["Stability"]], on=key, how="inner",
        suffixes=("_rg", "_conf"))
    if len(merged) < 0.99 * len(num):
        raise SystemExit("analytical and numerical full-plane grids do not "
                         "cover the same points; refusing to plot")

    mts, mhs, g = regular_grid(merged, ["Stability_rg", "Stability_conf", "S_exact"])
    st_rg = g["Stability_rg"]
    st_conf = g["Stability_conf"]

    # Difference categories (grid arrays for scatter).
    MM, HH = np.meshgrid(mts, mhs, indexing="ij")
    breakdown = ((st_rg == 3) & (g["S_exact"] <= 0)) | (st_rg == 4) | (st_conf == 4)
    rg_more_stable = (st_conf == 3) & (st_rg == 2)
    rg_less_stable = ((st_conf == 2) & (st_rg == 3)) | ((st_conf == 1) & (st_rg != 1))
    other = (st_rg != st_conf) & ~breakdown & ~rg_more_stable & ~rg_less_stable

    fig, ax = plt.subplots(figsize=(7.4, 6.4))
    ax.set_facecolor("0.94")  # agreement between the two classifiers
    for mask, color, label in (
            (rg_more_stable, "#2166ac", "RG-improved more stable"),
            (rg_less_stable, "#e08214", "RG-improved less stable"),
            (other, "0.35", "other disagreement"),
            (breakdown, "black", "breakdown / out of validity")):
        if mask.any():
            ax.plot(HH[mask], MM[mask], linestyle="none", marker="s",
                    markersize=2.4, color=color, label=label)
    ax.plot(BENCH_MH, BENCH_MT, marker="*", color="black", markersize=13,
            linestyle="none", markeredgecolor="white", markeredgewidth=0.5,
            zorder=6)

    ax.set_xlabel(r"$M_h$ [GeV]")
    ax.set_ylabel(r"$M_t$ [GeV]")
    ax.set_xlim(mhs[0], mhs[-1])
    ax.set_ylim(mts[0], mts[-1])
    ax.legend(loc="upper left", fontsize=7.5, framealpha=0.9, frameon=True,
              edgecolor="0.8", markerscale=3)

    n_diff = int((st_rg != st_conf).sum())
    print(f"differing points: {n_diff} of {st_rg.size}")
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
