# Example datasets

Pre-computed outputs of the stability classifiers, produced with the
production settings (`dt = 0.1`, `N = 2048`) by
`scripts/run_phase_diagram.py`. They back the figures in `figures/` and the
README; both can be regenerated with the commands below (about 4 minutes on
12 cores).

## `numerical_sm_region.csv` (12,705 rows)

RG-improved ansatz classification over the SM-relevant window
`Mt ∈ [155, 185] GeV`, `Mh ∈ [112, 138] GeV`, step 0.25 GeV:

```bash
python scripts/run_phase_diagram.py --mode numerical \
    --mt-min 155 --mt-max 185 --mh-min 112 --mh-max 138 --step 0.25 \
    --output data/numerical_sm_region.csv
```

Columns: `Mt, Mh_calc, Stability, S_exact, S_approx, S_kinetic,
S_potential, S_threshold, mu_inst, lambda_min` (schema in
[../docs/numerical-method.md](../docs/numerical-method.md)).

## `numerical_full_plane.csv` (63,001 rows)

RG-improved classification over the full overview plane
`Mt, Mh in [0, 250] GeV` at **1 GeV resolution**:

```bash
python scripts/run_phase_diagram.py --mode numerical     --mt-min 0 --mt-max 250 --mh-min 0 --mh-max 250 --step 1.0     --precision 8 --output data/numerical_full_plane.csv
```

(~10 min on 12 cores.) The 1 GeV resolution is a documented choice for the
global overview figure (features of interest are several GeV wide); the
original 0.25 GeV full-plane scan is 1,002,001 points and several hours —
the same command with `--step 0.25` reproduces it. The phenomenological
window is committed at 0.25 GeV in `numerical_sm_region.csv`.

## `analytical_full_plane.csv` (63,001 rows)

Strict conformal classification on the identical full-plane grid (same
command with `--mode analytical`). Used by the classification difference
map (figure 05).

## `analytical_sm_region.csv` (38,801 rows)

Strict conformal estimator over the wider window
`Mt ∈ [140, 200] GeV`, `Mh ∈ [105, 145] GeV`, step 0.25 GeV:

```bash
python scripts/run_phase_diagram.py --mode analytical \
    --mt-min 140 --mt-max 200 --mh-min 105 --mh-max 145 --step 0.25 \
    --output data/analytical_sm_region.csv
```

Columns: `Mt, Mh_calc, Stability, S_approx`.

## Status codes

| code | meaning |
|------|---------|
| 1 | stable (`lambda_eff` never negative up to the Planck scale) |
| 2 | metastable (action above the age-of-the-universe threshold) |
| 3 | unstable (action below the threshold; includes ansatz-breakdown points with `S_exact <= 0`) |
| 4 | perturbativity lost before the Planck scale |

Diagnostics (RG trajectories, action decomposition curves) for the figures
are generated on demand by the plotting scripts via
`build/dump_diagnostics` and cached under `results/` (not committed).
