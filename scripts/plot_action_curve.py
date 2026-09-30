#!/usr/bin/env python
"""Figure 08: the ansatz action S(R) and its kinetic/potential decomposition
at the benchmark point, with the golden-section minimum marked.

The horizontal axis R is the scale of the fixed Fubini-Lipatov profile
phi(r) = sqrt(2/|lambda_R|) * 2R/(R^2 + r^2), in units of the inverse
zero-crossing scale 1/mu1. The kinetic part is the closed-form conformal
result; the potential part (Simpson integral) carries the entire RG
improvement. The minimum of S(R) defines the reported action S*.

The physical meaning of R: the optimal radius corresponds to a physical
bubble radius R/mu1 (annotated in the figure), and the profile tip reaches
phi(0) = 2 sqrt(2/|lambda_R|) mu1 / R.
"""

import argparse
import os
import subprocess
import sys

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_style import apply_style

GEV_INV_TO_M = 1.973269804e-16  # hbar*c in GeV^-1 m
PLANCK_LENGTH_M = 1.616255e-35


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary",
                        default=os.path.join("build", "dump_diagnostics" + (".exe" if os.name == "nt" else "")))
    parser.add_argument("--data-dir", default="results/diagnostics")
    parser.add_argument("--output", default="figures/08_action_vs_R.png")
    args = parser.parse_args()

    apply_style()
    d = os.path.join(args.data_dir, "mh125.1_mt173.1")
    curve = os.path.join(d, "action_curve.csv")
    if not os.path.isfile(curve):
        os.makedirs(d, exist_ok=True)
        subprocess.run([args.binary, "--mh", "125.1", "--mt", "173.1",
                        "--out-dir", d], check=True, capture_output=True)
    df = pd.read_csv(curve)
    meta = pd.read_csv(os.path.join(d, "meta.csv")).set_index("key")["value"]
    mu1 = float(meta["mu1"])
    R_opt = float(meta["R_opt"])

    i_min = df["S_total"].idxmin()

    # Physical bubble radius of the optimal profile; the profile tip is
    # phi(0) = 2 sqrt(2/|lambda_R|) mu1 / R with lambda_R from the kinetic
    # term at the minimum (16 pi^2 / (3 |lambda_R|) = S_kinetic).
    import numpy as np
    R_phys_m = R_opt / mu1 * GEV_INV_TO_M
    R_in_lp = R_phys_m / PLANCK_LENGTH_M
    lam_R = 16.0 * np.pi ** 2 / (3.0 * df["S_kinetic"][i_min])
    phi0 = 2.0 * (2.0 / abs(lam_R)) ** 0.5 * mu1 / R_opt

    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    ax.plot(df["R_over_mu1"], df["S_kinetic"], color="#2166ac", lw=1.4,
            label=r"$S_{\rm kin} = 16\pi^2/(3|\lambda_R|)$ (closed form)")
    ax.plot(df["R_over_mu1"], -df["S_potential"], color="#ef8a62", lw=1.4,
            label=r"$|S_{\rm pot}|$ (Simpson integral)")
    ax.plot(df["R_over_mu1"], df["S_total"], color="black", lw=1.8,
            label=r"$S(R) = S_{\rm kin} + S_{\rm pot}$")

    s_star = df["S_total"][i_min]
    r_star = df["R_over_mu1"][i_min]
    ax.axvline(r_star, color="black", lw=0.8, linestyle=":", alpha=0.8)
    ax.annotate(rf"$S_\ast = {s_star:.1f}$" "\n"
                rf"$R_\ast/\mu_1 = {r_star:.1e}$",
                xy=(r_star, s_star), xytext=(r_star * 1.5, s_star * 1.05),
                fontsize=8.5, arrowprops=dict(arrowstyle="->", lw=0.8))

    ax.set_xscale("log")
    ax.set_xlabel(r"$R$  [units of $1/\mu_1$]")
    ax.set_ylabel("action contribution")
    ax.legend(frameon=False, fontsize=8, loc="upper center")
    ax.text(0.02, 0.02,
            rf"optimal profile: $R_{{\rm phys}} = R_\ast/\mu_1 \approx "
            rf"{R_phys_m:.1e}$ m $\approx {R_in_lp:.0f}\,\ell_{{\rm Pl}}$;  "
            rf"$\phi(0) \approx {phi0:.1e}$ GeV",
            transform=ax.transAxes, fontsize=7.5, va="bottom")

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
