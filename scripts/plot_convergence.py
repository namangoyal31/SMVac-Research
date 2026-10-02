#!/usr/bin/env python
"""Figure 10: numerical convergence of the production pipeline at the
benchmark point (125.1, 173.1) GeV.

Panel (a): the optimized ansatz action S* as a function of the RK4 step dt
(in t = ln(mu^2), log scale). The production step dt = 0.1 (marked)
reproduces the finest run to ~1e-6 relative (annotated).

Panel (b): |S(R*) - S(R*; N_max)| as a function of the number of Simpson
points N. The action is converged to the ~1e-9 relative round-off floor by
N ~ 256; production uses N = 2048. The dashed guide shows the O(N^-4)
Simpson scaling for orientation; below the floor the difference is
dominated by floating-point cancellation, not truncation.

Data generated on the fly by build/convergence_scan (not committed).
"""

import argparse
import os
import subprocess
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from io import StringIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_style import apply_style


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary",
                        default=os.path.join("build", "convergence_scan" + (".exe" if os.name == "nt" else "")))
    parser.add_argument("--output", default="figures/10_numerical_convergence.png")
    args = parser.parse_args()

    apply_style()
    dt_csv = subprocess.run([args.binary, "--mode", "dt"], check=True,
                            capture_output=True, text=True).stdout
    quad_csv = subprocess.run([args.binary, "--mode", "quad"], check=True,
                              capture_output=True, text=True).stdout
    dt_df = pd.read_csv(StringIO(dt_csv))
    quad_df = pd.read_csv(StringIO(quad_csv))

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10.4, 4.3))

    # (a) RK4 step convergence.
    from matplotlib.ticker import ScalarFormatter
    S_ref = dt_df["S_trial"].iloc[-1]
    ax_a.plot(dt_df["dt"], dt_df["S_trial"], marker="o", color="black",
              lw=1.0, markersize=4.5)
    ax_a.axvline(0.1, color="#c62828", lw=1.0, linestyle="--")
    ax_a.annotate(r"production $\Delta t = 0.1$", xy=(0.1, S_ref),
                  xytext=(0.11, S_ref + 0.012), fontsize=8, color="#c62828")
    d_prod = abs(dt_df.loc[dt_df["dt"] == 0.1, "S_trial"].iloc[0] - S_ref)
    ax_a.annotate(rf"$|S(0.1) - S(0.0125)|/S \approx {d_prod / S_ref:.0e}$",
                  xy=(0.03, S_ref - 0.02), fontsize=8)
    ax_a.set_xscale("log")
    ax_a.set_xticks([0.4, 0.2, 0.1, 0.05, 0.025, 0.0125])
    ax_a.set_xticklabels(["0.4", "0.2", "0.1", "0.05", "0.025", "0.0125"],
                         fontsize=8)
    yfmt = ScalarFormatter(useOffset=False)
    yfmt.set_powerlimits((-2, 2))
    ax_a.yaxis.set_major_formatter(yfmt)
    ax_a.set_xlabel(r"RK4 step $\Delta t$  [$t = \ln\mu^2$]")
    ax_a.set_ylabel(r"$S_\ast$ at $(125.1,\ 173.1)$ GeV")
    ax_a.set_title("(a) RG-step convergence", fontsize=9.5)
    pad = (dt_df["S_trial"].max() - dt_df["S_trial"].min()) * 4 + 1e-3
    ax_a.set_ylim(S_ref - pad, dt_df["S_trial"].max() + pad)

    # (b) quadrature convergence to the round-off floor.
    S_quad_ref = quad_df["S_at_Ropt"].iloc[-1]
    diff = np.abs(quad_df["S_at_Ropt"] - S_quad_ref)
    # Exclude the reference point itself (difference = 0 is not plottable
    # on a log axis).
    nz = diff > 0
    ax_b.plot(quad_df["N"][nz], diff[nz], marker="o", color="black", lw=1.0,
              markersize=4.5)
    ns = np.array([64.0, 16384.0])
    guide = diff.iloc[0] * (ns / ns[0]) ** (-4)
    ax_b.plot(ns, guide, linestyle="--", color="0.5", lw=0.9,
              label=r"$\mathcal{O}(N^{-4})$ guide")
    ax_b.axhline(1e-9 * S_quad_ref, color="#c62828", lw=0.9, linestyle=":")
    ax_b.annotate(r"round-off floor $\sim 10^{-9} S$",
                  xy=(200, 1.3e-9 * S_quad_ref), fontsize=8, color="#c62828")
    ax_b.axvline(2048, color="#1f77b4", lw=1.0, linestyle="--")
    ax_b.annotate(r"production $N = 2048$", xy=(2048, diff.max() * 0.4),
                  fontsize=8, color="#1f77b4")
    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.set_xlabel("Simpson points $N$")
    ax_b.set_ylabel(r"$|S(N) - S(N_{\rm max})|$ at $R_\ast$")
    ax_b.legend(frameon=False, fontsize=8)
    ax_b.set_title("(b) quadrature convergence", fontsize=9.5)

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
