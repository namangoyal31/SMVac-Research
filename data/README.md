# Datasets

Pre-computed outputs of the stability classifiers, produced with the
production settings (`dt = 0.1`, `N = 2048`) by
`scripts/run_phase_diagram.py`. **Every grid point in every file was
actually evaluated by the pipeline — no dataset is interpolated from a
coarser one**, and figure captions state the resolution they were built
from. The resolution-convergence between the grids is quantified by
`python scripts/check_resolution_convergence.py` (see
`docs/validation.md`).

| file | estimator | window | step | points | size |
|---|---|---|---|---|---|
| `numerical_full_plane_0p5GeV.csv` | RG-improved | Mh, Mt in [0, 250] GeV | 0.5 GeV | 251,001 | 15 MB |
| `analytical_full_plane_0p5GeV.csv` | conformal | Mh, Mt in [0, 250] GeV | 0.5 GeV | 251,001 | 4.9 MB |
| `numerical_zoom_0p1GeV.csv` | RG-improved | Mh [112, 138], Mt [155, 185] | 0.1 GeV | 78,561 | 5.5 MB |
| `analytical_zoom_0p1GeV.csv` | conformal | Mh [112, 138], Mt [155, 185] | 0.1 GeV | 78,561 | 1.7 MB |
| `numerical_boundary_0p05GeV.csv` | RG-improved | Mh [110, 140], Mt [162, 182] | 0.05 GeV | 241,001 | 18 MB |
| `numerical_zoom_0p1GeV_c6p0.csv` | RG-improved | Mh [112, 138], Mt [155, 185] | 0.1 GeV, c6 = 0 | 78,561 | 5.5 MB |
| `analytical_boundary_0p05GeV.csv` | conformal | Mh [110, 140], Mt [162, 182] | 0.05 GeV | 241,001 | 5.5 MB |
| `numerical_full_plane_1p0GeV.csv` | RG-improved | Mh, Mt in [0, 250] GeV | 1.0 GeV | 63,001 | 3.5 MB |
| `analytical_full_plane_1p0GeV.csv` | conformal | Mh, Mt in [0, 250] GeV | 1.0 GeV | 63,001 | 1.0 MB |
| `numerical_zoom_0p25GeV.csv` | RG-improved | Mh [112, 138], Mt [155, 185] | 0.25 GeV | 12,705 | 1.1 MB |
| `analytical_zoom_0p25GeV.csv` | conformal | Mh [112, 138], Mt [155, 185] | 0.25 GeV | 12,705 | 0.4 MB |

**Canonical production datasets**: the 0.5 GeV full plane (overview
figures), the 0.1 GeV zoom (phenomenological figures) and the 0.05 GeV
boundary strip (boundary/breakdown/difference diagnostics). The 1.0 GeV
full plane and 0.25 GeV zoom are retained as resolution-convergence
references - they are the coarse members of the convergence study in
`docs/validation.md`, not redundant copies.

**Boundary-window choice.** The 0.05 GeV strip covers
`Mh in [110, 140] x Mt in [162, 182] GeV` because the 0.25/0.1 GeV scans
show that *all* classification structure (stability boundary, breakdown
strip, threshold boundary over the full Mh range of the zoom) lies inside
that rectangle; outside it the classification is locally constant on
>= 1 GeV scales, where 0.05 GeV sampling would add cost without
information.

## Regeneration commands

```bash
# Full plane, 0.5 GeV (~28 min numerical / ~1 min analytical on 12 cores)
python scripts/run_phase_diagram.py --mode numerical     --mt-min 0 --mt-max 250 --mh-min 0 --mh-max 250 --step 0.5     --precision 8 --output data/numerical_full_plane_0p5GeV.csv
python scripts/run_phase_diagram.py --mode analytical     --mt-min 0 --mt-max 250 --mh-min 0 --mh-max 250 --step 0.5     --precision 8 --output data/analytical_full_plane_0p5GeV.csv

# Phenomenological zoom, 0.1 GeV (~12 min / ~0.5 min)
python scripts/run_phase_diagram.py --mode numerical     --mt-min 155 --mt-max 185 --mh-min 112 --mh-max 138 --step 0.1     --precision 8 --output data/numerical_zoom_0p1GeV.csv
python scripts/run_phase_diagram.py --mode analytical     --mt-min 155 --mt-max 185 --mh-min 112 --mh-max 138 --step 0.1     --precision 8 --output data/analytical_zoom_0p1GeV.csv

# Boundary strip, 0.05 GeV (~62 min / ~1 min)
python scripts/run_phase_diagram.py --mode numerical     --mt-min 162 --mt-max 182 --mh-min 110 --mh-max 140 --step 0.05     --precision 8 --output data/numerical_boundary_0p05GeV.csv
python scripts/run_phase_diagram.py --mode analytical     --mt-min 162 --mt-max 182 --mh-min 110 --mh-max 140 --step 0.05     --precision 8 --output data/analytical_boundary_0p05GeV.csv

# Coarser reference grids (same commands with --step 1.0 / 0.25)
```

Output precision is 8 significant digits (`--precision`), which preserves
the pipeline's own convergence level (~1e-6 relative) while bounding file
size.

## Columns and status codes

Numerical mode: `Mt, Mh_calc, Stability, S_trial, S_conformal, S_kinetic,
S_potential, S_threshold, mu_inst, lambda_min, method_flag, c6`;
analytical mode: `Mt, Mh_calc, Stability, S_conformal, method_flag`
(schema in [../docs/numerical-method.md](../docs/numerical-method.md)).
`method_flag`: 0 OK, 1 kinetic shortcut, 2 fence (no action computed,
`S_trial` = NaN), 3 ansatz failed (non-positive action), 4 perturbativity
lost. Status 0 always accompanies flags 2 and 3: those points are method
failures, not physical verdicts. The `c6` column records the
Planck-operator coefficient used (default 1).

| code | meaning |
|------|---------|
| 1 | stable (`lambda_eff` never negative up to the Planck scale) |
| 2 | metastable (action above the age-of-the-universe threshold) |
| 3 | unstable (action below the threshold; includes ansatz-breakdown points with `S_exact <= 0`) |
| 4 | perturbativity lost before the Planck scale |

Diagnostics (RG trajectories, action decomposition curves, convergence
traces) for the figures are generated on demand by
`build/dump_diagnostics` and `build/convergence_scan` and cached under
`results/` (not committed).
