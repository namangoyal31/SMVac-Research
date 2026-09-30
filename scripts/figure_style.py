"""Shared visual language for all SMVac-Research figures.

Physical colour semantics (matching the original SMVac plots and the
conventions of the vacuum-stability literature):

    stable          -> green
    metastable      -> yellow
    unstable        -> red
    non-perturbative / outside validity -> black

Continuous quantities (actions, fractional differences) use sequential or
diverging heatmaps instead; see docs/figures.md for the rationale.
"""

import matplotlib

matplotlib.use("Agg")

# Physical category colours.
COL_STABLE = "#3a9c3a"       # green
COL_METASTABLE = "#f5d327"   # yellow
COL_UNSTABLE = "#c62828"     # red
COL_INVALID = "#111111"      # black (non-perturbative / outside validity)
COL_BREAKDOWN = "#1a1a1a"    # ansatz-breakdown overlay markers (black dots)
COL_ACCENT = "#1f77b4"       # neutral accent for curves/markers

STATUS_COLORS = {
    1: COL_STABLE,
    2: COL_METASTABLE,
    3: COL_UNSTABLE,
    4: COL_INVALID,
}
STATUS_LABELS = {
    1: "stable",
    2: "metastable",
    3: "unstable",
    4: "non-perturbative",
}
BREAKDOWN_LABEL = r"ansatz breakdown ($S \leq 0$)"

# Benchmark point used throughout the code and reference tables.
BENCH_MH, BENCH_MT = 125.1, 173.1
# PDG 2022 central values and 1-sigma uncertainties (GeV).
PDG_MH, PDG_MH_ERR = 125.20, 0.11
PDG_MT, PDG_MT_ERR = 172.57, 0.29

RC = {
    "figure.dpi": 110,
    "savefig.dpi": 300,
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.linewidth": 0.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.frameon": False,
    "mathtext.fontset": "dejavusans",
    "axes.grid": False,
}


def apply_style():
    matplotlib.rcParams.update(RC)
