#!/usr/bin/env python
"""Figure: RG running of the effective quartic coupling.

Plots lambda_eff(mu) for three representative points of the (Mh, Mt) plane
(unstable, physical/metastable, and near-boundary metastable), showing the
zero crossing that drives vacuum decay and how its depth sets the action.

Trajectory data is generated with build/dump_diagnostics (cached under
results/diagnostics/), so the figure is always reproducible from the binary.
"""

import argparse
import os
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

POINTS = [
    ("115.0", "180.0", r"$(M_h,M_t)=(115,\ 180)$ GeV — unstable", "#b2182b"),
    ("125.1", "173.1", r"$(M_h,M_t)=(125.1,\ 173.1)$ GeV — metastable", "#2166ac"),
    ("134.75", "176.5", r"$(M_h,M_t)=(134.75,\ 176.5)$ GeV — near boundary", "#ef8a62"),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary",
                        default=os.path.join("build", "dump_diagnostics" + (".exe" if os.name == "nt" else "")))
    parser.add_argument("--data-dir", default="results/diagnostics")
    parser.add_argument("--output", default="figures/lambda_running.png")
    args = parser.parse_args()

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for mh, mt, label, color in POINTS:
        d = os.path.join(args.data_dir, f"mh{mh}_mt{mt}")
        traj = os.path.join(d, "trajectory.csv")
        if not os.path.isfile(traj):
            os.makedirs(d, exist_ok=True)
            subprocess.run([args.binary, "--mh", mh, "--mt", mt, "--out-dir", d],
                           check=True, capture_output=True)
        df = pd.read_csv(traj)
        ax.plot(df["scale_mu"], df["lambda_eff"], color=color, lw=1.6, label=label)

    ax.axhline(0.0, color="black", lw=0.8, alpha=0.7)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\mu$ [GeV]")
    ax.set_ylabel(r"$\lambda_{\mathrm{eff}}(\mu)$")
    ax.set_ylim(-0.06, 0.18)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.annotate(r"$\lambda_{\rm eff} < 0$ at high scales $\Rightarrow$ vacuum decay",
                xy=(1e7, -0.035), fontsize=8, color="black")

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
