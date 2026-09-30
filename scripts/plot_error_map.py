#!/usr/bin/env python
"""Figure: fractional difference between the RG-improved ansatz action and
the strict conformal estimate, (S_exact - S_approx)/S_approx, over the SM
region. The conformal estimate is exact when the couplings do not run over
the bubble; the map shows where the RG improvement is essential (near the
stability boundary) and where it is negligible (deep metastability).

Data source: data/numerical_sm_region.csv (produced by scripts/run_phase_diagram.py).
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_csv", nargs="?", default="data/numerical_sm_region.csv")
    parser.add_argument("--output", default="figures/error_map.png")
    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)
    valid = (df["S_exact"] > 0) & (df["S_approx"] > 0)
    rel = np.where(valid, (df["S_exact"] - df["S_approx"]) / df["S_approx"], np.nan)

    mts = np.sort(df["Mt"].unique())
    mhs = np.sort(df["Mh_calc"].unique())
    order = np.lexsort((df["Mh_calc"], df["Mt"]))
    Z = rel[order].reshape(len(mts), len(mhs))
    status = df["Stability"].to_numpy()[order].reshape(len(mts), len(mhs))

    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    vmax = np.nanpercentile(np.abs(Z), 98)
    mesh = ax.pcolormesh(mhs, mts, Z, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                         shading="auto", rasterized=True)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r"$(S_{\rm exact} - S_{\rm approx}) / S_{\rm approx}$")

    # The absolute stability boundary (status 1 region masked in the map)
    # for context.
    if (status == 1).any():
        ax.contour(mhs, mts, (status == 1).astype(float), levels=[0.5],
                   colors="black", linewidths=1.2)
        ax.plot([], [], color="black", lw=1.2, label=r"$\lambda_{\min} = 0$ boundary")
        ax.legend(frameon=False, fontsize=8, loc="upper left")

    ax.plot(125.1, 173.1, marker="*", color="black", markersize=12,
            markeredgecolor="white", markeredgewidth=0.5, linestyle="none")
    ax.annotate("SM point", xy=(125.1, 173.1), xytext=(126.5, 174.5), fontsize=8)

    ax.set_xlabel(r"$M_h$ [GeV]")
    ax.set_ylabel(r"$M_t$ [GeV]")
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
