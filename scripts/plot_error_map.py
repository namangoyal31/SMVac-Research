#!/usr/bin/env python
"""Figure 06: fractional difference between the RG-improved ansatz action and
the strict conformal estimate over the phenomenological window,

    delta(S) = (S_exact - S_approx) / S_approx,

evaluated at every point where both actions are defined (the metastable
region). The strict conformal estimate is exact when the couplings do not
run over the bubble; the map shows where the RG improvement materially
changes the action (near the absolute-stability boundary, where the profile
probes scales over which the couplings vary strongly) and where it is
negligible (deep in the metastable region). This is a continuous quantity,
so a diverging heatmap is used instead of the discrete green/yellow/red
phase colours.
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
    parser.add_argument("--output", default="figures/06_action_fractional_difference.png")
    args = parser.parse_args()

    apply_style()
    df = pd.read_csv("data/numerical_zoom_0p1GeV.csv")
    sel = ((df["Mt"] >= ZOOM_REGION["mt"][0]) & (df["Mt"] <= ZOOM_REGION["mt"][1]) &
           (df["Mh_calc"] >= ZOOM_REGION["mh"][0]) & (df["Mh_calc"] <= ZOOM_REGION["mh"][1]))
    mts, mhs, g = regular_grid(df[sel], ["S_exact", "S_approx", "Stability", "lambda_min"])

    valid = (g["S_exact"] > 0) & (g["S_approx"] > 0)
    Z = np.where(valid, (g["S_exact"] - g["S_approx"]) / g["S_approx"], np.nan)

    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    vmax = np.nanpercentile(np.abs(Z), 98)
    mesh = ax.pcolormesh(mhs, mts, Z, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                         shading="auto", rasterized=True)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r"$(S_{\rm exact} - S_{\rm approx}) / S_{\rm approx}$")

    # Absolute-stability boundary as the lambda_min = 0 level set of the
    # computed field (visualization-only interpolation between grid points).
    lam = np.where(g["Stability"] != 4, g["lambda_min"], np.nan)
    if np.isfinite(lam).any() and (lam > 0).any() and (lam < 0).any():
        ax.contour(mhs, mts, lam, levels=[0.0], colors="black", linewidths=1.2)
        ax.plot([], [], color="black", lw=1.2,
                label=r"$\lambda_{\min} = 0$ (absolutely stable)")
        ax.legend(frameon=False, fontsize=8, loc="upper left")

    ax.plot(BENCH_MH, BENCH_MT, marker="*", color="black", markersize=12,
            markeredgecolor="white", markeredgewidth=0.5, linestyle="none")
    ax.annotate("SM point", xy=(BENCH_MH, BENCH_MT),
                xytext=(BENCH_MH + 1.4, BENCH_MT + 1.4), fontsize=8)

    ax.set_xlabel(r"$M_h$ [GeV]")
    ax.set_ylabel(r"$M_t$ [GeV]")
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
