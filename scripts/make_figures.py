#!/usr/bin/env python
"""Regenerates every committed figure from the committed example data and the
compiled binaries. Run after building (cmake --build build):

    python scripts/make_figures.py

Figures (written to figures/):
    phase_diagram.png  - stability classification + action contours (SM region)
    lambda_running.png - lambda_eff(mu) for three representative points
    action_curve.png   - S(R) decomposition at the benchmark point
    error_map.png      - fractional difference S_exact vs S_approx
"""

import argparse
import subprocess
import sys

STEPS = [
    ("phase diagram", [sys.executable, "scripts/plot_phase_diagram.py",
                       "data/numerical_sm_region.csv",
                       "--output", "figures/phase_diagram.png", "--region", "sm"]),
    ("lambda running", [sys.executable, "scripts/plot_running.py"]),
    ("action curve", [sys.executable, "scripts/plot_action_curve.py"]),
    ("error map", [sys.executable, "scripts/plot_error_map.py",
                   "data/numerical_sm_region.csv",
                   "--output", "figures/error_map.png"]),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", choices=[s[0].replace(" ", "_") for s in STEPS])
    args = parser.parse_args()

    for name, cmd in STEPS:
        if args.only and name.replace(" ", "_") != args.only:
            continue
        print(f"== {name} ==")
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
