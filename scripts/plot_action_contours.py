#!/usr/bin/env python
"""Figure 04: magnitude of the RG-improved Fubini-Lipatov action over the
phenomenological window.

Filled map of log10 S_trial (continuous quantity -> heatmap, viridis scale),
with labeled constant-action contours. The action is only defined where the
running quartic coupling turns negative; in the absolutely stable region
(lower right) no bounce exists (hatched), and along the stability boundary
isolated points are blank where the fixed ansatz breaks down (S <= 0; see
figure 07 and docs/limitations.md).
"""

import argparse
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_style import apply_style, BENCH_MH, BENCH_MT
from plot_phase_diagram import regular_grid

ZOOM_REGION = dict(mh=(112.0, 138.0), mt=(155.0, 185.0))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", default="figures/04_action_contours.png")
    args = parser.parse_args()

    apply_style()
    df = pd.read_csv("data/numerical_zoom_0p1GeV.csv")
    sel = ((df["Mt"] >= ZOOM_REGION["mt"][0]) & (df["Mt"] <= ZOOM_REGION["mt"][1]) &
           (df["Mh_calc"] >= ZOOM_REGION["mh"][0]) & (df["Mh_calc"] <= ZOOM_REGION["mh"][1]))
    mts, mhs, g = regular_grid(df[sel], ["S_trial", "Stability"])

    # The action is only meaningful below the kinetic-shortcut/fence
    # region: rows with S_trial >= 4.9e5 are either the kinetic shortcut
    # (returned without the potential integral) or the search-bracket fence
    # (S = 1e100, no radius with lambda_R < 0). Both are classified
    # metastable, but the *value* is not a full action calculation; they are
    # blanked here and shown explicitly in figure 07.
    ACTION_MASK = 4.9e5
    valid = (np.isin(g["Stability"], [2, 3]) &
             (g["S_trial"] > 0) & (g["S_trial"] < ACTION_MASK))
    logS = np.where(valid, np.log10(np.where(g["S_trial"] > 0, g["S_trial"], 1.0)),
                    np.nan)

    fig, ax = plt.subplots(figsize=(7.0, 5.4))
    mesh = ax.pcolormesh(mhs, mts, logS, cmap="viridis", shading="auto",
                         rasterized=True)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r"$\log_{10} S_{\rm trial}$  ($S$ dimensionless)")

    levels = [450, 600, 800, 1200, 2000, 4000, 8000, 16000, 32000, 64000, 130000]
    cs = ax.contour(mhs, mts, np.where(valid, np.where(g["S_trial"] > 0,
                                                       g["S_trial"], np.nan), np.nan),
                    levels=levels, colors="k", linewidths=0.55, alpha=0.9)
    ax.clabel(cs, fmt=lambda v: f"S={v:.0f}", fontsize=6.3, inline=True,
              colors="white")

    # Blank strip adjacent to the stability boundary: kinetic-shortcut and
    # fence points (no full action integral; see figure 07).
    noblank = (np.isin(g["Stability"], [2, 3]) & (g["S_trial"] >= ACTION_MASK))
    if noblank.any():
        MM, HH = np.meshgrid(mts, mhs, indexing="ij")
        ax.plot(HH[noblank], MM[noblank], linestyle="none", marker="s",
                markersize=2.4, color="#fdae6b",
                label="shortcut/fence: no full action integral")
        ax.legend(loc="upper right", fontsize=7.5, framealpha=0.9,
                  frameon=True, edgecolor="0.8")

    # Absolutely stable region: no bounce exists.
    stable = g["Stability"] == 1
    if stable.any():
        ax.contourf(mhs, mts, stable.astype(float), levels=[0.5, 1.5],
                    colors="white", hatches=["////"])
        ax.contour(mhs, mts, stable.astype(float), levels=[0.5],
                   colors="black", linewidths=1.2)
        ax.text(0.97, 0.05, r"$\lambda_{\min} \geq 0$:" "\n" "absolutely stable, no bounce",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5)

    ax.plot(BENCH_MH, BENCH_MT, marker="*", color="black", markersize=13,
            linestyle="none", markeredgecolor="white", markeredgewidth=0.5,
            zorder=6)
    ax.annotate("SM point", xy=(BENCH_MH, BENCH_MT),
                xytext=(BENCH_MH - 3.4, BENCH_MT + 1.1), fontsize=8)

    ax.set_xlabel(r"$M_h$ [GeV]")
    ax.set_ylabel(r"$M_t$ [GeV]")

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
