#!/usr/bin/env python
"""Figure 06: fractional difference between the trial-profile action and the
strict conformal estimate over the phenomenological window,

    delta(S) = (S_trial - S_conformal) / S_conformal,

at 0.1 GeV sampling, shown for two values of the Planck-suppressed
operator coefficient: c6 = 0 (left) isolates the effect of RG improvement
inside the trial family; c6 = 1 (right, the production default) adds the
assumed Planck operator, which dominates near the boundary. A diverging
red-blue heatmap is used: this is a continuous quantity, not a phase
classification.
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


def frac_panel(ax, df, label):
    mts, mhs, g = regular_grid(df, ["S_trial", "S_conformal", "Stability"])
    valid = (g["S_trial"] > 0) & (g["S_conformal"] > 0)
    Z = np.where(valid, (g["S_trial"] - g["S_conformal"]) / g["S_conformal"], np.nan)

    vmax = np.nanpercentile(np.abs(Z), 98)
    mesh = ax.pcolormesh(mhs, mts, Z, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                         shading="auto", rasterized=True)
    ax.set_xlabel(r"$M_h$ [GeV]")
    ax.set_ylabel(r"$M_t$ [GeV]")
    ax.set_title(label, fontsize=9.5)
    return mesh


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--c6-0-csv", default="data/numerical_zoom_0p1GeV_c6p0.csv")
    parser.add_argument("--c6-1-csv", default="data/numerical_zoom_0p1GeV.csv")
    parser.add_argument("--output", default="figures/06_action_fractional_difference.png")
    args = parser.parse_args()

    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.7), sharey=True)

    panels = []
    for ax, path, label in (
            (axes[0], args.c6_0_csv, r"$c_6 = 0$ (RG improvement only)"),
            (axes[1], args.c6_1_csv, r"$c_6 = 1$ (production default)")):
        df = pd.read_csv(path)
        sel = ((df["Mt"] >= ZOOM_REGION["mt"][0]) & (df["Mt"] <= ZOOM_REGION["mt"][1]) &
               (df["Mh_calc"] >= ZOOM_REGION["mh"][0]) & (df["Mh_calc"] <= ZOOM_REGION["mh"][1]))
        panels.append(frac_panel(ax, df[sel], label))

    # Level-set stability boundary on both panels.
    for ax, path in ((axes[0], args.c6_0_csv), (axes[1], args.c6_1_csv)):
        df = pd.read_csv(path)
        sel = ((df["Mt"] >= ZOOM_REGION["mt"][0]) & (df["Mt"] <= ZOOM_REGION["mt"][1]) &
               (df["Mh_calc"] >= ZOOM_REGION["mh"][0]) & (df["Mh_calc"] <= ZOOM_REGION["mh"][1]))
        mts, mhs, g = regular_grid(df[sel], ["Stability", "lambda_min"])
        lam = np.where(g["Stability"] != 4, g["lambda_min"], np.nan)
        if np.isfinite(lam).any() and (lam > 0).any() and (lam < 0).any():
            ax.contour(mhs, mts, lam, levels=[0.0], colors="black", linewidths=1.1)
    axes[0].plot([], [], color="black", lw=1.1,
                 label=r"$\lambda_{\min} = 0$ (absolutely stable)")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    axes[0].plot(BENCH_MH, BENCH_MT, marker="*", color="black", markersize=11,
                 markeredgecolor="white", markeredgewidth=0.5, linestyle="none")
    axes[1].plot(BENCH_MH, BENCH_MT, marker="*", color="black", markersize=11,
                 markeredgecolor="white", markeredgewidth=0.5, linestyle="none")

    vmin = min(p_.get_clim()[0] for p_ in panels)
    vmax = max(p_.get_clim()[1] for p_ in panels)
    for p_ in panels:
        p_.set_clim(vmin, vmax)
    cbar = fig.colorbar(panels[1], ax=axes, pad=0.015)
    cbar.set_label(r"$(S_{\rm trial} - S_{\rm conformal}) / S_{\rm conformal}$")

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
