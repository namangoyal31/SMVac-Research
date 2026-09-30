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

## [1.0.0] - 2026-07-01 (original SMVac repository)

- Initial public release: 3-loop SM RGE running, RG-improved 1-loop
  potential, Fubini-Lipatov and RG-improved conformal-profile bounce
  estimators, phase-diagram pipeline.
