#!/usr/bin/env python
"""Phase-diagram figures: stability classification of the SM vacuum in the
(Mh, Mt) plane, from the RG-improved Fubini-Lipatov ansatz.

Physical colour semantics (as in the original SMVac plots):
    green  = stable        (lambda_eff(mu) >= 0 up to the Planck scale)
    yellow = metastable    (lambda_eff < 0, action above the cosmic threshold)
    red    = unstable      (action below the threshold)
    black  = non-perturbative / outside validity (status 4)

Points where the fixed conformal ansatz breaks down (classified unstable
with a non-positive action, adjacent to the stability boundary) are overlaid
as black dots: they are a documented artifact of the method, not physical
instabilities (docs/limitations.md, item 12).

Regions:
    --region full : the whole scanned plane 0-250 GeV x 0-250 GeV
                    (data/numerical_full_plane.csv, 1 GeV resolution)
    --region zoom : the phenomenological window Mh 112-138, Mt 155-185 GeV
                    (data/numerical_sm_region.csv, 0.25 GeV resolution)
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
from matplotlib.patches import Ellipse, Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_style import (apply_style, STATUS_COLORS, STATUS_LABELS,
                          BENCH_MH, BENCH_MT, PDG_MH, PDG_MT, PDG_MH_ERR,
                          PDG_MT_ERR)

ZOOM_REGION = dict(mh=(112.0, 138.0), mt=(155.0, 185.0))


def regular_grid(df, cols):
    """Reshape the Mt-major scan CSV into (Mt, Mh) grids for each column."""
    mts = np.sort(df["Mt"].unique())
    mhs = np.sort(df["Mh_calc"].unique())
    shape = (len(mts), len(mhs))
    order = np.lexsort((df["Mh_calc"], df["Mt"]))
    grids = {}
    for col in cols:
        vals = df[col].to_numpy(float)[order]
        if vals.size != shape[0] * shape[1]:
            raise SystemExit(
                f"CSV is not a complete regular grid: {vals.size} values for "
                f"a {shape[0]}x{shape[1]} grid ({shape[0]*shape[1]} expected). "
                "Refuse to plot silently-substituted data.")
        grids[col] = vals.reshape(shape)
    return mts, mhs, grids


def classification_map(ax, mts, mhs, g, contour_levels, ellipses, mt_step):
    """Draw one classification map on the given axes."""
    status = g["Stability"]
    cmap = ListedColormap([STATUS_COLORS[i] for i in (1, 2, 3, 4)])
    norm = BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)
    mh_edges = np.append(mhs - mt_step / 2, mhs[-1] + mt_step / 2)
    mt_edges = np.append(mts - mt_step / 2, mts[-1] + mt_step / 2)
    ax.pcolormesh(mh_edges, mt_edges, status, cmap=cmap, norm=norm,
                  shading="flat", rasterized=True)

    # Constant-action contours of the RG-improved ansatz action (where the
    # action is defined).
    S = np.where(g["S_exact"] > 0, g["S_exact"], np.nan)
    if contour_levels:
        cs = ax.contour(mhs, mts, S, levels=contour_levels, colors="k",
                        linewidths=0.6, alpha=0.65)
        ax.clabel(cs, fmt=lambda v: f"S={v:.0f}", fontsize=6.5, inline=True)

    # Boundary of the genuine (non-artifact) unstable region.
    unb = np.where((g["Stability"] == 3) & (g["S_exact"] > 0), 1.0, 0.0)
    if 0 < unb.sum() < unb.size:
        ax.contour(mhs, mts, unb, levels=[0.5], colors="k",
                   linewidths=0.9, linestyles="dashed")

    # Ansatz-breakdown points (unstable with S <= 0): black dots.
    bd = (g["Stability"] == 3) & (g["S_exact"] <= 0)
    if bd.any():
        MM, HH = np.meshgrid(mts, mhs, indexing="ij")
        ax.plot(HH[bd], MM[bd], linestyle="none", marker="o", markersize=2.2,
                markerfacecolor="black", markeredgecolor="black", zorder=5)

    # Benchmark point; PDG ellipses where they resolve at this scale.
    ax.plot(BENCH_MH, BENCH_MT, marker="*", color="black", markersize=13,
            linestyle="none", zorder=6, markeredgecolor="white",
            markeredgewidth=0.5)
    if ellipses:
        for factor, style in ((1, "--"), (2, ":")):
            ax.add_patch(Ellipse((PDG_MH, PDG_MT),
                                 width=2 * factor * PDG_MH_ERR,
                                 height=2 * factor * PDG_MT_ERR, fill=False,
                                 edgecolor="black", linestyle=style,
                                 linewidth=1.0, zorder=5))

    ax.set_xlabel(r"$M_h$ [GeV]")
    ax.set_ylabel(r"$M_t$ [GeV]")
    ax.set_xlim(mhs[0], mhs[-1])
    ax.set_ylim(mts[0], mts[-1])

    handles = [Patch(facecolor=STATUS_COLORS[s], label=STATUS_LABELS[s])
               for s in (1, 2, 3, 4)]
    if bd.any():
        handles.append(Line2D([], [], linestyle="none", marker="o", markersize=4,
                              markerfacecolor="black", markeredgecolor="black",
                              label="ansatz breakdown ($S\\leq 0$)"))
    handles.append(Line2D([], [], linestyle="none", marker="*", markersize=10,
                          markerfacecolor="black", markeredgecolor="white",
                          label="benchmark (125.1, 173.1) GeV"))
    if ellipses:
        handles.append(Line2D([], [], linestyle="--", color="black",
                              label="PDG 2022 $1\\sigma$"))
    ax.legend(handles=handles, loc="upper left", fontsize=7.5, framealpha=0.9,
              frameon=True, edgecolor="0.8")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--region", choices=["full", "zoom"], required=True)
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    apply_style()
    if args.region == "full":
        df = pd.read_csv("data/numerical_full_plane.csv")
        mts, mhs, g = regular_grid(df, ["Stability", "S_exact"])
        mt_step = mts[1] - mts[0]
        output = args.output or "figures/01_full_phase_diagram.png"
        levels = [450, 600, 800, 1200, 2000, 4000, 10000]
        fig, ax = plt.subplots(figsize=(7.4, 6.4))
        classification_map(ax, mts, mhs, g, levels, ellipses=False,
                           mt_step=mt_step)
    else:
        df = pd.read_csv("data/numerical_sm_region.csv")
        sel = ((df["Mt"] >= ZOOM_REGION["mt"][0]) & (df["Mt"] <= ZOOM_REGION["mt"][1]) &
               (df["Mh_calc"] >= ZOOM_REGION["mh"][0]) & (df["Mh_calc"] <= ZOOM_REGION["mh"][1]))
        mts, mhs, g = regular_grid(df[sel], ["Stability", "S_exact"])
        mt_step = mts[1] - mts[0]
        output = args.output or "figures/02_phenomenological_zoom.png"
        levels = [450, 600, 800, 1200, 2000, 4000, 8000, 16000, 32000]
        fig, ax = plt.subplots(figsize=(7.0, 5.6))
        classification_map(ax, mts, mhs, g, levels, ellipses=True,
                           mt_step=mt_step)

    fig.tight_layout()
    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    fig.savefig(output, dpi=300)
    print(f"Saved {output}")


if __name__ == "__main__":
    main()
