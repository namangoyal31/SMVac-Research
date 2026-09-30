#!/usr/bin/env python
"""Figure 09: RG running at the benchmark point.

Panel (a): the effective quartic coupling lambda_eff(mu) for three
representative points of the (Mh, Mt) plane. The zero crossing and its
depth control the vacuum classification; the unstable point (low Mh, high
Mt) crosses early and deeply, the physical point crosses at
mu1 ~ 6e10 GeV, and the near-boundary point barely turns negative.

Panel (b): all SM couplings actually evolved by the code at the benchmark
point (g1 GUT-normalized, g2, g3, y_t, lambda). The top Yukawa drives
lambda down; the running stops at the Planck scale where the couplings are
frozen (docs/theory.md).
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
from figure_style import apply_style

POINTS = [
    ("115.0", "180.0", r"$(115,\ 180)$ GeV (unstable)", "#c62828"),
    ("125.1", "173.1", r"$(125.1,\ 173.1)$ GeV (metastable)", "#1f77b4"),
    ("134.75", "176.5", r"$(134.75,\ 176.5)$ GeV (near boundary)", "#e68a2e"),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary",
                        default=os.path.join("build", "dump_diagnostics" + (".exe" if os.name == "nt" else "")))
    parser.add_argument("--data-dir", default="results/diagnostics")
    parser.add_argument("--output", default="figures/09_rg_running.png")
    args = parser.parse_args()

    apply_style()
    trajs = {}
    for mh, mt, _, _ in POINTS:
        d = os.path.join(args.data_dir, f"mh{mh}_mt{mt}")
        traj = os.path.join(d, "trajectory.csv")
        if not os.path.isfile(traj):
            os.makedirs(d, exist_ok=True)
            subprocess.run([args.binary, "--mh", mh, "--mt", mt,
                            "--out-dir", d], check=True, capture_output=True)
        trajs[(mh, mt)] = pd.read_csv(traj)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10.4, 4.4))

    # (a) lambda_eff for three points.
    for (mh, mt, label, color) in POINTS:
        df = trajs[(mh, mt)]
        ax_a.plot(df["scale_mu"], df["lambda_eff"], color=color, lw=1.5,
                  label=label)
    ax_a.axhline(0.0, color="black", lw=0.8, alpha=0.7)
    ax_a.set_xscale("log")
    ax_a.set_xlabel(r"$\mu$ [GeV]")
    ax_a.set_ylabel(r"$\lambda_{\rm eff}(\mu)$")
    ax_a.set_ylim(-0.06, 0.18)
    ax_a.legend(frameon=False, fontsize=7.5, loc="upper right", handlelength=1.6)
    ax_a.set_title("(a) effective quartic coupling", fontsize=9.5)

    # (b) all evolved couplings at the benchmark point.
    df = trajs[("125.1", "173.1")]
    ax_b.plot(df["scale_mu"], df["yt"], color="#762a83", lw=1.5,
              label=r"$y_t$")
    ax_b.plot(df["scale_mu"], df["g3"], color="#1b7837", lw=1.5,
              label=r"$g_3$")
    ax_b.plot(df["scale_mu"], df["g2"], color="#2166ac", lw=1.5,
              label=r"$g_2$")
    ax_b.plot(df["scale_mu"], df["g1"], color="#e08214", lw=1.5,
              label=r"$g_1$ (GUT norm.)")
    ax_b.plot(df["scale_mu"], df["lambda_raw"], color="black", lw=1.5,
              label=r"$\lambda$")
    ax_b.axhline(0.0, color="black", lw=0.6, alpha=0.5)
    ax_b.set_xscale("log")
    ax_b.set_xlabel(r"$\mu$ [GeV]")
    ax_b.set_ylabel("coupling")
    ax_b.legend(frameon=False, fontsize=7.5, loc="upper left", handlelength=1.6)
    ax_b.set_title(r"(b) SM couplings at $(125.1,\ 173.1)$ GeV", fontsize=9.5)

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
