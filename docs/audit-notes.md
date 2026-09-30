# Code Audit Notes

*Audit performed before the v2.0 reorganization of this repository. The
corresponding source snapshot and the full working analysis notes remain in the
public [`SMVac`](https://github.com/namangoyal31/SMVac) repository; this
document records the findings and their disposition in the present codebase.*

## Scope

Every source file, application, test, script, configuration file, reference
dataset, and document in the repository was reviewed, with the physics and
numerics of the core library checked line by line against the intended
calculation. The baseline was established by building the library with
GCC 6.3.0 (MinGW-w64, `-std=c++17 -O2`) and running the full test suite; all
frozen reference values reproduce to `~1e-15` relative error on that
configuration.

## What the code actually computes

The library estimates the electroweak vacuum decay action of the Standard
Model as follows:

1. **Matching.** NNLO electroweak matching conditions (linearized around the
   central masses) fix the SM couplings at the top-mass scale.
2. **Running.** The 3-loop SM beta functions for `(g1, g2, g3, yt, yb, ytau,
   lambda)` are integrated with fixed-step RK4 (`dt = 0.1` in the evolution
   variable) from `Mt` up to the Planck scale.
3. **Potential.** The RG-improved 1-loop effective potential is evaluated as
   `V(phi) = lambda_eff(phi) * phi^4 / 4`, with couplings run to `mu = phi`
   and the Coleman–Weinberg correction from top, W and Z loops.
4. **Action.** The bounce action is evaluated on the Fubini–Lipatov conformal
   profile with running-coupling amplitude, over bubble radii `R` optimized by
   golden-section search; the potential integral is done with composite
   Simpson quadrature.
5. **Classification.** Points are labelled stable / metastable / unstable by
   comparing the optimized action with the age-of-the-universe threshold.

The correct characterization of the method is a **semi-analytical bounce
action estimator based on the Fubini–Lipatov profile ansatz with numerical RG
running and numerical potential integration**. It is *not* a solution of the
Euclidean bounce equation: no shooting, relaxation, or collocation is
performed, and the profile family is constrained to the conformal form with a
single optimized scale. Because the true bounce minimizes the action over all
profiles, the reported action is a variational **upper bound** on the exact
bounce action. Earlier draft wording ("exact numerical bounce solver") was
incorrect and has been removed everywhere.

## Findings and disposition

| # | Finding | Severity | Disposition |
|---|---------|----------|-------------|
| 1 | The legacy monolith `apps/solver_numerical.cpp` duplicated the entire library, and its `--analytical` mode only changed the output *filename* while still computing the numerical classifier — analytical-labelled data from that binary would have been mislabelled | high | Removed; the library (whose analytical mode is implemented correctly) is the single source of truth |
| 2 | `RGEHelper::get_params` linearly interpolated **only** `lambda` between RG table points; all other couplings were held at the left bracket point (staircase evaluation of the potential) | medium | Fixed: all couplings interpolated; the shift of the benchmark actions is quantified in `docs/validation.md` |
| 3 | The RG evolution variable is `t = ln(mu^2)` and all beta functions are implemented as `dX/d ln(mu^2)` (i.e. half the standard `dX/d ln(mu)` rates). This is self-consistent throughout but nonstandard and easy to misread as a factor-2 bug | documentation | Documented in `docs/numerical-method.md` and in code comments |
| 4 | The potential is the high-field form `V = lambda_eff(phi) phi^4/4`: the tree-level mass term is dropped and the false vacuum is `phi = 0`, not the electroweak vacuum at 246 GeV. Intentional (conformal-ansatz) approximation | documentation | Documented in `docs/theory.md`; listed under known limitations |
| 5 | A Planck-suppressed `phi^6` operator with hard-coded coefficient `c6 = 1` regularizes the action near the Planck scale. No literature provenance was established for this specific coefficient | medium (theory) | Promoted to a named constant with documentation; a sensitivity study is listed under known limitations |
| 6 | The metastability threshold hid the universe age in the magic constant `1.179e44` (which equals `t_U * v` with `t_U = 10 Gyr` in GeV^-1) | medium | Refactored to named, derived constants; the 10 Gyr choice is documented |
| 7 | The decay-rate prefactor uses the lambda zero-crossing scale; the `1/R` dependence of the prefactor is not included in the action minimization (cf. the complete-lifetime treatment of Andreassen, Frost & Schwartz 2018) | known approximation | Documented in `docs/theory.md` and under known limitations |
| 8 | Action evaluation returns the kinetic term alone when `16 pi^2 / (3 |lambda_R|) > 5e5`, skipping the potential integral | low | Documented; affects only `|lambda| <~ 3e-4`, where the action exceeds the threshold by orders of magnitude regardless |
| 9 | Dead code: `V1_loop`, `rk4_adaptive_step` (defined but never called; documentation claimed adaptive stepping was used), constant `SIGMA_PLAN` | low | Removed; the fixed-step RK4 convergence is quantified in `docs/validation.md` |
| 10 | `apps/generate_phase_diagram.cpp` never checked `ofstream::is_open()`; with a missing output directory all writes failed silently while the program printed "Done" | medium | Fixed: explicit error on failure; directories created by the workflow scripts |
| 11 | Documentation overclaims relative to the code: "adaptive RK4" (fixed `dt = 0.1` is used), "OpenMP parallelization" (no OpenMP pragmas exist; parallelism is Python-level multiprocessing), "exact numerical bounce solver" | high (presentation) | Corrected throughout |
| 12 | The existing tests compare against the code's own frozen output (regression tests), and the RK4 convergence study in the working notes used a toy ODE (`dy/dt = -ky`), not the actual RGE system | high (validation) | Real validation added: analytic limits, genuine convergence studies of the actual pipeline, and literature checks — see `docs/validation.md` |
| 13 | The `g3` beta function contains a 4-loop term with numeric coefficient `2472.28` whose provenance could not be established from the code or the cited literature | flagged | Kept (removing it would silently change the physics); its impact is measured: benchmark actions shift by ~0.10% when the term is dropped. Listed under known limitations |
| 14 | `dAlphas = (alpha3_at_Mz - 0.1184)/0.0007` is identically zero: alpha_s is hard-coded to its central value | low | Documented; the parameterization is retained for future scans |
| 15 | The repository contained debris from the legacy research repository: documentation of an unrelated third-party package, working notes with machine-specific paths, duplicated Python runners, an empty `CMakeLists.txt`, an empty `LICENSE`, and a `CITATION.cff` with a placeholder URL | medium | Removed; provenance of the migrated code retained in `docs/CODE_PROVENANCE.md` |
| 16 | The NNLO matching conditions are a linearized fit around `(Mh, Mt) = (125.15, 173.34) GeV`; the `Mh`-dependence of `delta lambda` is not implemented | medium | Documented; the matching is tested against published central values in `tests/unit`, and the approximation is listed under known limitations |
| 17 | Points with negative bounce action (ansatz breakdown far outside the conformal regime) are classified as unstable | known behaviour | Documented; the kinetic/potential decomposition is now exported so such points can be diagnosed |

## Verified claims carried forward

Two quantitative statements in the previous README were re-verified by direct
computation and are retained:

* Deep metastability, `(Mh, Mt) = (115, 180) GeV`: the RG-improved action
  agrees with the conformal estimate `8 pi^2 / (3 |lambda_min|)` to
  `-0.015%` (claimed: `-0.01%`).
* Near the boundary, `(Mh, Mt) = (134.75, 176.5) GeV`: the two differ by
  `+32.1%` (claimed: `32.1%`), confirming that the running-coupling
  corrections are essential near the metastability boundary.
