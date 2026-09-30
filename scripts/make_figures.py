#!/usr/bin/env python
"""Regenerates every committed figure from the committed datasets and the
compiled binaries. Run after building (cmake --build build):

    python scripts/make_figures.py

Figures (written to figures/, documented in docs/figures.md):
    01_full_phase_diagram            full 0-250 GeV classification map
    02_phenomenological_zoom         SM-window classification at 0.25 GeV
    03_experimental_point            benchmark point vs stability boundary
    04_action_contours               action magnitude over the SM window
    05_conformal_vs_rg_difference    classification difference map
    06_action_fractional_difference  (S_exact - S_approx)/S_approx heatmap
    07_ansatz_breakdown              validity/breakdown map (full plane)
    08_action_vs_R                   S(R) decomposition at the SM point
    09_rg_running                    RG running of the couplings
    10_numerical_convergence         dt and N convergence of the pipeline

Datasets: data/numerical_full_plane.csv and data/numerical_sm_region.csv
(committed; see data/README.md). Diagnostics (trajectories, action curves,
convergence traces) are generated on the fly by build/dump_diagnostics and
build/convergence_scan and cached under results/ (gitignored).
"""

import argparse
import subprocess
import sys

STEPS = [
    ("01_full_phase_diagram", [sys.executable, "scripts/plot_phase_diagram.py", "--region", "full"]),
    ("02_phenomenological_zoom", [sys.executable, "scripts/plot_phase_diagram.py", "--region", "zoom"]),
    ("03_experimental_point", [sys.executable, "scripts/plot_experimental_point.py"]),
    ("04_action_contours", [sys.executable, "scripts/plot_action_contours.py"]),
    ("05_conformal_vs_rg_difference", [sys.executable, "scripts/plot_conformal_vs_rg.py"]),
    ("06_action_fractional_difference", [sys.executable, "scripts/plot_error_map.py"]),
    ("07_ansatz_breakdown", [sys.executable, "scripts/plot_validity.py"]),
    ("08_action_vs_R", [sys.executable, "scripts/plot_action_curve.py"]),
    ("09_rg_running", [sys.executable, "scripts/plot_running.py"]),
    ("10_numerical_convergence", [sys.executable, "scripts/plot_convergence.py"]),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", choices=[name for name, _ in STEPS])
    args = parser.parse_args()

    for name, cmd in STEPS:
        if args.only and name != args.only:
            continue
        print(f"== {name} ==")
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
