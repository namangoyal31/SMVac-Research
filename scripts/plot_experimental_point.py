#!/usr/bin/env python
"""Figure 03: the benchmark point (Mh, Mt) = (125.1, 173.1) GeV relative to
the calculated stability structure.

A focused view of the phenomenological window: phase classification from
the RG-improved ansatz (green/yellow/red), the absolute-stability boundary,
the benchmark point, and its quantitative distance to the stability
boundary. Annotations (S*, threshold, zero-crossing scale) are read from
the diagnostics produced by build/dump_diagnostics.
"""

import argparse
import os
import subprocess
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_style import (apply_style, BENCH_MH, BENCH_MT, PDG_MH, PDG_MT,
                          PDG_MH_ERR, PDG_MT_ERR)
from plot_phase_diagram import classification_map, regular_grid

WINDOW = dict(mh=(121.0, 133.0), mt=(167.0, 179.0))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary",
                        default=os.path.join("build", "dump_diagnostics" + (".exe" if os.name == "nt" else "")))
    parser.add_argument("--data-dir", default="results/diagnostics")
    parser.add_argument("--output", default="figures/03_experimental_point.png")
    args = parser.parse_args()

    apply_style()
    df = pd.read_csv("data/numerical_boundary_0p05GeV.csv")
    sel = ((df["Mt"] >= WINDOW["mt"][0]) & (df["Mt"] <= WINDOW["mt"][1]) &
           (df["Mh_calc"] >= WINDOW["mh"][0]) & (df["Mh_calc"] <= WINDOW["mh"][1]))
    mts, mhs, g = regular_grid(df[sel], ["Stability", "S_exact", "lambda_min", "S_threshold"])
    mt_step = mts[1] - mts[0]

    # Diagnostics at the benchmark point (regenerate via dump_diagnostics).
    d = os.path.join(args.data_dir, "mh125.1_mt173.1")
    meta_path = os.path.join(d, "meta.csv")
    if not os.path.isfile(meta_path):
        os.makedirs(d, exist_ok=True)
        subprocess.run([args.binary, "--mh", "125.1", "--mt", "173.1",
                        "--out-dir", d], check=True, capture_output=True)
    meta = pd.read_csv(meta_path).set_index("key")["value"]
    S_star = float(meta["S_star"])
    S_th = float(meta["S_threshold"])
    mu1 = float(meta["mu1"])

    # Distance to the absolute-stability boundary along the benchmark row.
    row = np.argmin(np.abs(mts - BENCH_MT))
    stable_cols = np.where(g["Stability"][row] == 1)[0]
    mh_crit = mhs[stable_cols[0]] if len(stable_cols) else np.nan
    d_stab = mh_crit - BENCH_MH

    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    classification_map(ax, mts, mhs, g, contour_levels=[600, 1200, 4000],
                       ellipses=True, mt_step=mt_step)

    # Quantitative annotation of the benchmark point.
    mu1_exp = int(np.floor(np.log10(mu1)))
    mu1_mant = mu1 / 10.0 ** mu1_exp
    ax.annotate(
        f"metastable:  $S_\\ast = {S_star:.0f} > S_{{\\rm th}} = {S_th:.0f}$\n"
        f"$\\lambda_{{\\rm eff}}$ crosses zero at $\\mu_1 = {mu1_mant:.1f}"
        f"\\times10^{{{mu1_exp}}}$ GeV"
        f"\n$\\Delta M_h$ to absolute stability $\\approx {d_stab:.1f}$ GeV",
        xy=(BENCH_MH, BENCH_MT), xytext=(BENCH_MH + 1.2, BENCH_MT + 3.2),
        fontsize=8.5, arrowprops=dict(arrowstyle="->", lw=0.9),
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.6", alpha=0.9))

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
