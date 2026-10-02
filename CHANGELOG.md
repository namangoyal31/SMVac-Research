# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2026-09-30

Upgraded research codebase, reorganized from the public
[namangoyal31/SMVac](https://github.com/namangoyal31/SMVac) repository after
a full audit (see `docs/audit-notes.md`).

### Fixed (physics)

- **Electroweak matching corrected to the cited NNLO references**: the
  previous `g2(Mt) = 0.65355` and `g1(Mt) = 0.46266` could not be traced to
  any published source and contradict Buttazzo et al. (2013), eqs. (58)-(59)
  (`g2(Mt) = 0.64779`, `g1(Mt) = sqrt(5/3)*0.35830 = 0.462563`). The
  correction lowers the benchmark actions by 2% (deep metastability) to 24%
  (near the stability boundary); the absolute-stability boundary moves to
  `Mh_crit = 129.2 GeV` at `Mt = 173.1` (published: ~128.6 GeV).
- All couplings are now interpolated in the RG lookup table (previously
  only `lambda`; the remaining couplings were staircase-held between table
  points).
- Frozen reference tables regenerated accordingly (`reference/v1.1/`, via
  the committed `apps/write_reference_tables.cpp`).

### Added

- Real validation suite (11 tests via CTest): analytic conformal limit
  (closed-form ansatz action reproduced to 1e-12..1e-16), RK4/Simpson
  convergence orders, production-pipeline convergence studies, matching vs
  published values, metastability-threshold derivation, literature
  comparison of the stability boundary. See `docs/validation.md`.
- CMake build system (library, five apps, CTest integration) and a GitHub
  Actions workflow running the full suite.
- Exported action decomposition (`S_kinetic`, `S_potential`), threshold,
  and diagnostic scales in `StabilityResult` and the scan CSVs.
- Reproducible workflow: `scripts/run_phase_diagram.py` (argparse-driven,
  failure-checked) and committed example datasets under `data/`.
- Publication figures with regeneration scripts
  (`scripts/make_figures.py`), including a dedicated category for the
  documented ansatz-breakdown points (`S <= 0`).
- Documentation set: `docs/theory.md`, `docs/numerical-method.md`,
  `docs/validation.md`, `docs/limitations.md`, `docs/references.md`,
  `docs/audit-notes.md`.

### Changed

- Metastability-threshold constants replaced by named, derived constants
  (the magic number `1.179e44` encoded `t_U = 10 Gyr` in GeV^-1 multiplied
  by the vev; threshold values change negligibly, ~2e-3).
- Output precision of scan CSVs raised to 12 significant digits.
- Bulk generated CSVs (65 MB) replaced by regenerable example datasets
  (2 MB); legacy working notes and unrelated files remain in the original
  repository only.

### Removed

- Duplicated legacy monolith solver (`apps/solver_numerical.cpp`) whose
  `--analytical` mode computed numerical results under an analytical
  filename.
- Dead code (`V1_loop`, the unused adaptive RK4 stepper, unused constants).
- Legacy constant naming cleaned up (`planck_mass`; see docs/audit-notes.md).

## [2.2.0] - 2026-10-02 (referee-correction pass)

Substantive scientific corrections implementing an independent referee
report (`SMVac_referee_report.md`, auditing HEAD a1bca57). All datasets and
figures regenerated; every grid point re-evaluated.

### Fixed (physics)
- **Z Coleman-Weinberg normalization (BUG)**: the Z mass combination used
  the GUT-normalized g1^2 instead of (3/5) g1^2 + g2^2
  (src/EffectivePotential.cpp). Post-fix benchmark values: SM-point trial
  action 2103.10 -> 2168.69 (+3.1%), mu1 6.22e10 -> 7.60e10 GeV,
  Mh_crit(173.1) 129.17 -> 129.05 GeV.
- **4-loop g3 sign**: 2472.28 traced to the pure-QCD beta_3(nf=6)
  (van Ritbergen, Vermaseren & Larin 1997); the code's own convention
  (-beta_i g^(2i+4), with -7, -26, +65/2 at nf=6) requires -2472.28. The
  previous +2472.28 had the wrong sign.

### Changed (claims)
- **"Variational upper bound on the true bounce action" retracted
  everywhere**: the bounce is a saddle point of the action, not a minimum
  over general field configurations; the code's quantity is now called the
  trial-profile action and no bound relation to the exact bounce action is
  claimed (docs/theory.md §5.1, with the Coleman-Glaser-Martin 1977
  scope statement).
- **"+30% RG improvement near the boundary" retracted**: with c6 = 0 the
  trial action matches the conformal estimate to <1% at every tested
  point; the near-boundary enhancement is the assumed c6 = 1 Planck
  operator. c6 is now a runtime parameter (--c6) and the sensitivity scan
  is committed (RESULTS.md, tests/physics/test_method_flags.cpp).
- **Literature comparison re-anchored**: the unverified "128.6 GeV
  (Buttazzo Fig. 1)" removed; tests/physics/test_boundary.cpp now uses
  Degrassi et al. (2012) eq. (2), 129.4 +- 1.0 GeV at the same inputs
  (Mt = 173.1, alpha_s = 0.1184), with like-for-like caveats stated.
- Renames: S_exact -> S_trial, S_approx -> S_conformal,
  classify_buttazzo -> classify_conformal (names that implied more than
  the calculation establishes).

### Added
- **method_flag column** in all CSVs and in StabilityResult: OK /
  KINETIC_SHORTCUT / FENCE / ANSATZ_FAILED / PERTURBATIVITY_LOST. Fence
  points get S_trial = NaN (no more 1e100 sentinels in output); points
  with non-positive trial action are status 0 (undetermined), never
  counted as metastable/unstable.
- tests/unit/test_lambda_eff_normalization.cpp (Z normalization,
  discriminating against the old formula), tests/physics/test_method_flags.cpp
  (flag taxonomy + c6 separation); suite extended 11 -> 13 tests.
- RESULTS.md: the defensible results after the corrections.

### Fixed (documentation)
- "16-decade" bracket corrected (+-8 in ln R = ~7 decades in R); M_W
  "implemented" claim corrected (it is not an input of the matching
  function); (115, 180) "deep metastability" label corrected (the code
  classifies it unstable); AFS arXiv number corrected (1707.08123 ->
  1707.08124); g1(Mt) value unified at 0.462563; gauge-dependence of mu1
  and the amplitude/potential mismatch of the trial family documented as
  open limitations; full-plane figure labeled an extended computational
  visualization (linearized matching).

## [2.1.0] - 2026-09-30

### Changed
- Phase-diagram resolution increased, with **all datasets actually
  recomputed** (no interpolation): full plane 1.0 -> **0.5 GeV**
  (251,001 points, ~28 min), phenomenological zoom 0.25 -> **0.1 GeV**
  (78,561 points, ~12 min), and a new dedicated boundary strip at
  **0.05 GeV** (`Mh ∈ [110, 140] × Mt ∈ [162, 182]`, 241,001 points,
  ~62 min). Datasets renamed to carry their resolution
  (`*_full_plane_0p5GeV.csv`, `*_zoom_0p1GeV.csv`,
  `*_boundary_0p05GeV.csv`); the coarse grids are kept as
  resolution-convergence references.
- Stability boundaries in the figures are now drawn as level sets of the
  per-point computed fields (`λ_min = 0`; `S = S_threshold`) instead of
  pixel edges of the classification; the interpolation between grid
  points is documented as purely visual.
- All derived figures (01-07) regenerated from the new-resolution data; a
  new figure 11 shows the breakdown and classification-difference strips
  at 0.05 GeV (both ~1 grid cell wide, ≲ 0.1 GeV).

### Added
- `scripts/check_resolution_convergence.py`: quantifies reproducibility
  (bit-identical common points), boundary-location convergence (the
  `λ_min = 0` crossing at Mt = 173 GeV is 128.965 GeV at every
  resolution), classification area fractions, disagreement fractions, and
  breakdown-count scaling. Documented in `docs/validation.md` §6.

## [2.0.1] - 2026-09-30

### Fixed
- `tests/unit/test_numerics`: the golden-section location check demanded
  1e-8 absolute precision — below the sqrt(eps*|f|) floating-point
  resolution of value-comparison minimization for the tested quadratic
  (~4e-8), which made the outcome platform-dependent (CI failed at
  2.1e-8). The test now verifies the algorithm's actual contracts (exact
  bracket contraction, function-value optimality, location within a
  documented 25x margin of the fp floor) plus a well-conditioned kink
  function localized to ~1e-14.

### Changed
- **License removed at the author's discretion**: the MIT LICENSE file and
  all license references are gone; the repository is copyright (c) 2026
  Naman Goyal, all rights reserved (see NOTICE). A copyright notice is not
  a license.
- Figure suite rebuilt (10 figures with the original green/yellow/red/black
  phase semantics, full-plane 0-250 GeV overview restored at 1 GeV
  resolution, difference and validity maps, convergence figure); see
  docs/figures.md. Full-plane datasets committed at 1 GeV resolution
  (63,001 points, `--precision 8`), generated by the documented commands.

## [1.0.0] - 2026-07-01 (original SMVac repository)

- Initial public release: 3-loop SM RGE running, RG-improved 1-loop
  potential, Fubini-Lipatov and RG-improved conformal-profile bounce
  estimators, phase-diagram pipeline.
