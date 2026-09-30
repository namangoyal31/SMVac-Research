#!/usr/bin/env python
"""Figure: ansatz action S(R) and its kinetic/potential decomposition at the
benchmark point, with the golden-section minimum marked. The kinetic term is
the closed-form conformal result; the potential integral carries the entire
RG improvement. In the pure-quartic limit S_kinetic/|S_potential| = 2.

Data source: build/dump_diagnostics (called automatically if missing).
"""

import argparse
import os
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary",
                        default=os.path.join("build", "dump_diagnostics" + (".exe" if os.name == "nt" else "")))
    parser.add_argument("--data-dir", default="results/diagnostics")
    parser.add_argument("--output", default="figures/action_curve.png")
    args = parser.parse_args()

    d = os.path.join(args.data_dir, "mh125.1_mt173.1")
    curve = os.path.join(d, "action_curve.csv")
    if not os.path.isfile(curve):
        os.makedirs(d, exist_ok=True)
        subprocess.run([args.binary, "--mh", "125.1", "--mt", "173.1", "--out-dir", d],
                       check=True, capture_output=True)
    df = pd.read_csv(curve)

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.plot(df["R_over_mu1"], df["S_kinetic"], color="#2166ac", lw=1.4,
            label=r"$S_{\rm kin} = 16\pi^2/(3|\lambda_R|)$ (closed form)")
    ax.plot(df["R_over_mu1"], -df["S_potential"], color="#ef8a62", lw=1.4,
            label=r"$|S_{\rm pot}|$ (Simpson integral)")
    ax.plot(df["R_over_mu1"], df["S_total"], color="black", lw=1.8,
            label=r"$S(R) = S_{\rm kin} + S_{\rm pot}$")

    i_min = df["S_total"].idxmin()
    s_star = df["S_total"][i_min]
    r_star = df["R_over_mu1"][i_min]
    ax.axvline(r_star, color="black", lw=0.8, linestyle=":", alpha=0.8)
    ax.annotate(rf"$S_\ast = {s_star:.1f}$", xy=(r_star, s_star),
                xytext=(r_star * 1.6, s_star * 1.02), fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=0.8))

    ax.set_xscale("log")
    ax.set_xlabel(r"$R$ [in units of $1/\mu_1$]")
    ax.set_ylabel("action contribution")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
