#!/usr/bin/env python
"""Resolution-convergence check for the phase-diagram datasets.

Compares every committed scan resolution and reports:

1. Reproducibility: S_exact at grid points common to all datasets must be
   bit-identical (the pipeline is deterministic; any difference would mean
   a dataset was produced by different code or settings).
2. Stability-boundary location at fixed Mt: grid-level last transition and
   the linearly interpolated lambda_min = 0 crossing (the interpolation is
   a diagnostic of the sampled field, not an extra calculation).
3. Classification area fractions (counts scale with grid density; area
   fractions are the convergent quantity).
4. Conformal-vs-RG classification disagreement (count and fraction).
5. Ansatz-breakdown point counts.
6. Classification of the experimental point (nearest grid point).

Run after (re)generating the datasets:

    python scripts/check_resolution_convergence.py
"""

import glob
import os
import sys

import numpy as np
import pandas as pd

# Common test points present on every grid (1.0, 0.5, 0.25, 0.1, 0.05 all
# contain integer mass values in their windows).
TEST_POINTS = [(125.0, 173.0), (120.0, 170.0), (130.0, 176.0)]
BENCH = (125.1, 173.1)


def load(pattern):
    files = glob.glob(pattern)
    return files[0] if files else None


def boundary_Mh_at_Mt(df, mt_value):
    """Grid-level and interpolated stability-boundary Mh at the given Mt."""
    row = df.iloc[np.argmin(np.abs(df["Mt"].to_numpy() - mt_value))]
    mt_actual = row["Mt"]
    d = df[np.isclose(df["Mt"], mt_actual)].sort_values("Mh_calc")
    # Grid-level boundary: last non-stable, non-perturbative point (status 4
    # regions at high Mh are not part of the stability boundary).
    nonstable = d[(d["Stability"] != 1) & (d["Stability"] != 4)]
    if nonstable.empty or (d["Stability"] == 1).all():
        return mt_actual, None, None
    mh_grid = nonstable["Mh_calc"].max()
    # lambda_min = 0 crossing by linear interpolation between the two
    # bracketing grid points (visualization/diagnostic only).
    lam = d["lambda_min"].to_numpy()
    mh = d["Mh_calc"].to_numpy()
    crossing = None
    for i in range(len(lam) - 1):
        if np.isfinite(lam[i]) and np.isfinite(lam[i + 1]) and lam[i] < 0 <= lam[i + 1]:
            t = -lam[i] / (lam[i + 1] - lam[i])
            crossing = mh[i] + t * (mh[i + 1] - mh[i])
            break
    return mt_actual, mh_grid, crossing


def status_fractions(df):
    n = len(df)
    return {s: df["Stability"].value_counts().get(s, 0) / n for s in (1, 2, 3, 4)}


def main():
    datasets = {}
    for name, pattern in [
        ("full 1.0 GeV", "data/numerical_full_plane_1p0GeV.csv"),
        ("full 0.5 GeV", "data/numerical_full_plane_0p5GeV.csv"),
        ("zoom 0.25 GeV", "data/numerical_zoom_0p25GeV.csv"),
        ("zoom 0.1 GeV", "data/numerical_zoom_0p1GeV.csv"),
        ("boundary 0.05 GeV", "data/numerical_boundary_0p05GeV.csv"),
    ]:
        f = load(pattern)
        if f:
            datasets[name] = pd.read_csv(f)

    print("# Resolution-convergence check\n")

    print("## 1. Reproducibility at common grid points\n")
    print("| point | " + " | ".join(datasets) + " |")
    print("|---|" + "---|" * len(datasets))
    for (mh, mt) in TEST_POINTS:
        vals = []
        for name, df in datasets.items():
            sel = df[(np.isclose(df["Mt"], mt)) & (np.isclose(df["Mh_calc"], mh))]
            vals.append(f"{sel['S_exact'].iloc[0]:.4f}" if len(sel) else "n/a")
        print(f"| ({mh}, {mt}) | " + " | ".join(vals) + " |")
    print()

    print("## 2. Absolute-stability boundary at Mt = 173 GeV\n")
    print("| dataset | grid-level last non-stable Mh | lambda_min = 0 crossing |")
    print("|---|---|---|")
    for name, df in datasets.items():
        mt_actual, mh_grid, crossing = boundary_Mh_at_Mt(df, 173.0)
        cg = f"{mh_grid:.2f}" if mh_grid is not None else "n/a"
        cc = f"{crossing:.3f}" if crossing is not None else "n/a"
        print(f"| {name} (Mt={mt_actual:.2f}) | {cg} | {cc} |")
    print()

    print("## 3. Classification area fractions\n")
    print("| dataset | stable | metastable | unstable | non-pert | breakdown pts |")
    print("|---|---|---|---|---|---|")
    for name, df in datasets.items():
        fr = status_fractions(df)
        nbd = int(((df["Stability"] == 3) & (df["S_exact"] <= 0)).sum())
        print(f"| {name} | {fr[1]:.4f} | {fr[2]:.4f} | {fr[3]:.4f} | "
              f"{fr[4]:.4f} | {nbd} |")
    print()

    print("## 4. Conformal-vs-RG classification disagreement\n")
    for num_pat, ana_pat, label in [
        ("data/numerical_full_plane_1p0GeV.csv", "data/analytical_full_plane_1p0GeV.csv", "full 1.0 GeV"),
        ("data/numerical_full_plane_0p5GeV.csv", "data/analytical_full_plane_0p5GeV.csv", "full 0.5 GeV"),
        ("data/numerical_zoom_0p1GeV.csv", "data/analytical_zoom_0p1GeV.csv", "zoom 0.1 GeV"),
        ("data/numerical_boundary_0p05GeV.csv", "data/analytical_boundary_0p05GeV.csv", "boundary 0.05 GeV"),
    ]:
        fn, fa = load(num_pat), load(ana_pat)
        if not (fn and fa):
            continue
        num = pd.read_csv(fn)
        ana = pd.read_csv(fa)
        key = ["Mt", "Mh_calc"]
        m = num[key + ["Stability"]].merge(ana[key + ["Stability"]], on=key,
                                           suffixes=("_rg", "_conf"))
        nd = int((m["Stability_rg"] != m["Stability_conf"]).sum())
        print(f"- {label}: {nd} of {len(m)} points differ ({nd/len(m):.4%})")
    print()

    print("## 5. Experimental point (nearest grid point)\n")
    print("| dataset | nearest (Mh, Mt) | status |")
    print("|---|---|---|")
    for name, df in datasets.items():
        idx = (np.abs(df["Mh_calc"] - BENCH[0]) +
               np.abs(df["Mt"] - BENCH[1])).idxmin()
        row = df.loc[idx]
        print(f"| {name} | ({row['Mh_calc']:.2f}, {row['Mt']:.2f}) | "
              f"{int(row['Stability'])} |")
    print()


if __name__ == "__main__":
    sys.exit(main())
