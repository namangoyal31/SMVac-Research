# Numerical Method

This document maps every physics equation onto the algorithm and the
implementation that evaluates it, states the numerical settings and their
validation status, and defines the output formats. The physics itself is
documented in [theory.md](theory.md).

---

## 1. Pipeline overview

```
 (Mh, Mt)
    |
    v
 [matching]            get_nnlo_matching            src/RGE.cpp
    |   MS-bar couplings (g1,g2,g3,yt,yb,ytau,lambda) at mu = Mt
    v
 [RG running]          rk4_single_step              src/RGE.cpp
    |   fixed-step RK4, dt = 0.1 in t = ln(mu^2), Mt -> M_Pl
    |   fills an RGEHelper table (linear interpolation in t)
    v
 [stability check]     get_lambda_eff               src/EffectivePotential.cpp
    |   lambda_eff(mu) = lambda(mu) + Delta_lambda_CW(mu)
    |   lambda_min, zero-crossing scale mu1
    |
    +--> lambda_min >= 0 ............................. status 1 (stable)
    +--> couplings non-perturbative .................. status 4
    |                                                   (flag PERTURBATIVITY_LOST)
    v
 [trial action]        evaluate_action_components   src/CanonicalBounce.cpp
    |   S(R) = 16 pi^2/(3|lambda_R|)  +  Simpson integral of 2 pi^2 r^3 V
    |   on the Fubini-Lipatov profile phi = sqrt(2/|lambda_R|) 2R/(R^2+r^2)
    v
 [radius optimization] find_minimum_action          src/CanonicalBounce.cpp
    |   golden-section search on log R, bracket +-8 in ln R (~7 decades
    |   in R), tol 1e-13
    v
 [classification]      classify_stability           src/CanonicalBounce.cpp
    |   no valid R ................. status 0 (flag FENCE,       S_trial = NaN)
    |   S* <= 0 .................... status 0 (flag ANSATZ_FAILED)
    |   S* > 4 ln(mu1 * t_U) ....... status 2 (metastable)
    |   else ....................... status 3 (unstable)
    v
 (StabilityResult / CSV)
```

The conformal estimator `classify_conformal` (src/FubiniLipatov.cpp) shares
the matching and running stages and stops after the stability check,
reporting the strict conformal action $S_{\rm conformal} =
8\pi^2/(3|\lambda_{\min}|)$.

## 2. Physics-to-code mapping

| Physics | Formulation | Algorithm | Implementation | Output |
|---|---|---|---|---|
| NNLO matching at $\mu = M_t$ (Buttazzo 2013, eqs. 55–60) | linearized fit; $M_h$-dependence via $M_h^2/2v^2$ | direct evaluation | `get_nnlo_matching` | initial `StandardModelParameters` |
| 3-loop SM RGEs, $\mathrm{d}X/\mathrm{d}\ln\mu^2$ | ODE system, no explicit $t$-dependence | classical RK4, fixed step | `rk4_single_step` + `RGEHelper::add_point` | RG trajectory table |
| RG-improved potential $V = \lambda_{\rm eff}\phi^4/4$ | $\mu = \phi$; CW correction with $\phi$-dependent masses; $m_Z^2 \supset g_Y^2 = (3/5) g_1^2$ | table lookup + linear interpolation | `V_eff`, `get_lambda_eff`, `RGEHelper::get_params` | potential values |
| Stability boundary | $\lambda_{\min} \ge 0$ | scan along trajectory | inside `classify_stability` / `classify_conformal` | status 1 vs 2/3 |
| Fubini–Lipatov profile $\phi_R = A\,2R/(R^2{+}r^2)$, $A^2 = 2/\|\lambda_R\|$ | fixed profile, scale $R$ | closed form | `integrand_u`, `get_pure_sm_lambda` | profile values |
| Kinetic action $S_{\rm kin} = 16\pi^2/(3\|\lambda_R\|)$ | analytic for the fixed profile | closed form | `evaluate_action_components` | `S_kinetic` |
| Potential action $S_{\rm pot} = 2\pi^2\int r^3 V\,\mathrm{d}r$ | compactified $u = \tanh x$, $r = Re^{\alpha x}$, $\alpha = 4$ | composite Simpson, $N = 2048$ | `integrand_u` + `numerics::simpson_integrate` | `S_potential` |
| Planck regularization $V \supset (c_6/6)(\phi/M_{\rm Pl})^2\phi^4$ | extra term in the integrand | — | `c6` parameter (default `kDefaultC6 = 1`; `--c6`) | (included in `S_potential`) |
| Action minimum $S_\ast$ | 1-parameter minimization in $\ln R$ | golden-section search, tol $10^{-13}$, bracket $\pm 8$ in $\ln R$ | `find_minimum_action` + `numerics::golden_section_search` | `S_trial` |
| Conformal estimate $8\pi^2/(3\|\lambda_{\min}\|)$ | closed form | — | `classify_stability`, `classify_conformal` | `S_conformal` |
| Metastability threshold $4\ln(\mu_1 t_U)$ | Coleman rate argument | closed form | `classify_stability` (constants in `RGE.hpp`) | `S_threshold`, status 2/3 |
| Method flag | sentinel/shortcut/failure bookkeeping | — | `classify_stability` (`MethodFlag` enum) | `method_flag` |

## 3. Numerical settings

| Setting | Value | Where | Validation |
|---|---|---|---|
| RK4 step $\Delta t$ (in $t = \ln\mu^2$) | 0.1 | `classify_stability`, `classify_conformal` | step-halving study: $\lambda_{\min}$ and $S_\ast$ converge to $\le 10^{-6}$ relative; see `tests/physics/test_convergence.cpp` |
| Quadrature points $N$ (Simpson, on $u \in (-1,1)$) | 2048 | `kDefaultQuadraturePoints` | at the $\sim 10^{-9}$ relative round-off floor for all $N$ tested (order $\sim N^{-4}$ verified on smooth integrands in `tests/unit/test_numerics.cpp`); the integrable $u \to -1$ edge singularity is validated against the analytic conformal limit to $\le 2\times10^{-12}$ |
| Golden-section tolerance on $\ln R$ | $10^{-13}$ | `kLogRadiusTolerance` | minimization error far below machine-relevant level; bracket $[\ln R_{\rm opt} \pm 8]$ (~7 decades in $R$) |
| Interpolation of the RG table | linear in $t$ | `RGEHelper::get_params` | error bounded by the table step; covered by the $\Delta t$ study |
| Running range | $\mu \in [M_t, M_{\rm Pl}]$ (production), extended to $[1\ \mathrm{GeV}, M_{\rm Pl}]$ in reference tables | — | couplings frozen at the range edges |
| Planck cap | $h \le M_{\rm Pl}$ in the integrand | `get_pure_sm_lambda` | documented approximation |
| Kinetic shortcut | returns $S_{\rm kin}$ alone if $S_{\rm kin} > 5\times10^5$ ($\|\lambda_R\| \lesssim 3.3\times10^{-4}$) | `kKineticShortcutThreshold` | flagged `KINETIC_SHORTCUT`; classification robust (value exceeds threshold by orders of magnitude) |

All frozen regression settings are reproduced by `tests/regression/`.

## 4. Output formats

### 4.1 `StabilityResult` (C++ struct, single point)

Defined in `include/SMVacuumDecay/CanonicalBounce.hpp`; fields are documented
in [theory.md §7](theory.md#7-interpretation-of-outputs).

### 4.2 Chunk CSV from `generate_phase_diagram`

Numerical mode:

```
Mt,Mh_calc,Stability,S_trial,S_conformal,S_kinetic,S_potential,S_threshold,mu_inst,lambda_min,method_flag,c6
```

Analytical mode:

```
Mt,Mh_calc,Stability,S_conformal,method_flag
```

`Mt` is the top-mass grid value, `Mh_calc` echoes the Higgs-mass grid value.
`Stability` uses the codes of [theory.md §6](theory.md#6-metastability-criterion):
1 stable, 2 metastable, 3 unstable, 4 perturbativity lost, **0 undetermined**
(always paired with `method_flag` FENCE or ANSATZ_FAILED — these points are
method failures, not physical verdicts; fence rows carry `S_trial = NaN`,
ansatz-failed rows keep their computed non-positive value). `method_flag`:
0 OK, 1 kinetic shortcut, 2 fence, 3 ansatz failed, 4 perturbativity lost.
`c6` records the Planck-operator coefficient used.

### 4.3 Reference tables (`reference/v1.1/`, generated by `apps/write_reference_tables.cpp`)

* `rge_reference.csv` — trajectory `scale_mu,g1,g2,g3,yt,lambda` for the
  benchmark point, from 1 GeV to $M_{\rm Pl}$.
* `potential_reference.csv` — `phi,lambda_eff,V` along the same trajectory.
* `benchmark_results*.json` — single-point actions (`Conformal_action`,
  `Trial_action`), statuses and method flags.

## 5. Reproducing a result

From a fresh clone (see the README for prerequisites):

```bash
cmake -S . -B build && cmake --build build
./build/benchmark_point 125.1 173.1               # single point, both estimators
python scripts/run_phase_diagram.py --mode numerical --mt-min 160 --mt-max 185 \
       --mh-min 110 --mh-max 140 --step 0.5       # SM-region scan
python scripts/plot_phase_diagram.py data/numerical_zoom_0p1GeV.csv \
       --output figures/02_phenomenological_zoom.png --region zoom  # visualize
ctest --test-dir build
```

The committed datasets and the figures in `figures/` are produced by
exactly these commands (see `scripts/make_figures.py`, `data/README.md`
and the README).
