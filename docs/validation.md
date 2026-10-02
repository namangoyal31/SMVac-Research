# Validation

This document records what has actually been validated, how, and with what
result. It distinguishes **analytic validation** (comparison with closed-form
results), **convergence studies** (discretization errors of the real
pipeline), **literature comparisons**, and **regression tests** (agreement
with frozen historical output). Nothing here is claimed without an
automated test that enforces it; the tests live in `tests/unit`,
`tests/physics`, and `tests/regression`.

Build and run everything with:

```bash
cmake -S . -B build && cmake --build build
ctest --test-dir build --output-on-failure
```

---

## 1. Analytic validation

### 1.1 Constant-coupling conformal limit (`tests/physics/test_fubini_limit.cpp`)

With a frozen coupling table ($\lambda = -0.01$ const., negligible Yukawa/gauge
content), the Fubini–Lipatov profile *is* the exact bounce and the action has
the closed form

$$S(R) = \frac{16\pi^2}{3|\lambda_R|} - \frac{8\pi^2\,|\lambda_{\rm eff}|}{3\lambda_R^2},$$

where $\lambda_R$ is the raw coupling used for the profile amplitude and
$\lambda_{\rm eff}$ the coupling inside the potential integrand. This is a
complete end-to-end check of the profile parametrization, the compactified
radial integral ($r = Re^{\alpha x}$, $u = \tanh x$), the Simpson quadrature,
and the closed-form kinetic term.

**Result:** the identity holds to $3\times10^{-16}$, $2\times10^{-14}$, and
$2\times10^{-12}$ relative error at profile scales $R = 10^{-3}, 10^{-2},
10^{-1}$ (in $1/\mu_1$ units). The strict Fubini–Lipatov value
$S = 8\pi^2/(3|\lambda_{\rm eff}|)$ is reproduced within the
$\lambda_{\rm eff}\ne\lambda_R$ amplitude correction. Scale invariance
($R$-independence of $S$) holds to $2\times10^{-12}$.

### 1.2 Matching conditions vs published values (`tests/unit/test_matching.cpp`)

The implemented matching, evaluated at the central masses
$(M_h, M_t) = (125.15, 173.34)$ GeV, compared with Buttazzo et al. (2013):

| coupling | implemented | published | difference |
|---|---|---|---|
| $\lambda(M_t)$ | 0.1260470 | 0.12604 | $7.0\times10^{-6}$ |
| $y_t(M_t)$ | 0.9369128 | 0.93690 | $1.3\times10^{-5}$ |
| $g_2(M_t)$ | 0.64779 | 0.64779 | 0 (set to the published value; see audit notes) |
| $g_3(M_t)$ | 1.1666 | 1.1666 | 0 |
| $g_1(M_t)$ | 0.462563 | $\sqrt{5/3}\times0.35830 = 0.462563$ | 0 |

### 1.3 Metastability threshold (`tests/unit/test_threshold.cpp`)

The age-of-the-universe constant is re-derived independently in the test
(10 Gyr × days × $\hbar^{-1}$) and agrees with the header constant to
$7\times10^{-4}$ (the residual is the 365- vs 365.25-day year). The threshold
$S_{\rm th} = 4\ln(\mu_1 t_U)$ at $\mu_1 = 10^{10}$ GeV equals 475.99,
inside the physically expected range $[458, 504]$ for $\mu_1 \in
[10^8, 10^{13}]$ GeV.

## 2. Convergence studies (`tests/physics/test_convergence.cpp`)

Studies of the **production pipeline** at the SM benchmark point
$(125.1, 173.1)$ GeV:

* **RGE step size.** $\Delta t \in \{0.2, 0.1, 0.05, 0.025\}$ gives
  $S_\ast = 2168.6962,\ 2168.6896,\ 2168.6879,\ 2168.6874$: successive
  differences shrink by a factor $\simeq 4$ per halving, and the production
  step $\Delta t = 0.1$ agrees with the four-times-finer run to
  $1.0\times10^{-6}$ relative. The RG evolution is therefore not a
  significant source of error at production settings.
* **Quadrature.** At the optimal radius, $S(R_\ast)$ for $N \in \{256,
  1024, 4096, 16384\}$ Simpson points agrees to $10^{-9}$ relative
  between $N = 2048$ and $N = 16384$. The differences between successive $N$
  are non-monotone at the $\sim10^{-9}$ relative level: the production
  setting $N = 2048$ already sits at the round-off floor (the potential
  integral involves large cancellations), so no meaningful order estimate
  is possible there. The theoretical $\mathcal{O}(N^{-4})$ Simpson
  convergence is verified separately on smooth integrands in
  `tests/unit/test_numerics.cpp` (observed order 4.00).
* **Radius minimization.** The golden-section search agrees with a 50,001-point
  grid scan of the same bracket ($\pm 8$ in $\ln R$, ~7 decades in $R$) to $1.2\times10^{-10}$ relative.

## 3. Literature comparison (`tests/physics/test_boundary.cpp`)

* **Physical point.** $(M_h, M_t) = (125.1, 173.1)$ GeV is classified
  **metastable** (status 2), matching the consensus that the SM vacuum
  lifetime exceeds the age of the universe.
* **Absolute stability boundary.** Bisection at $M_t = 173.1$ GeV locates
  $\lambda_{\min} = 0$ at $M_h^{\rm crit} = 129.05$ GeV. Anchor (verified
  against the published text): Degrassi et al. (2012), introduction
  eq. (2): $M_h > 129.4 + 1.4\frac{M_t - 173.1}{0.7} -
  0.5\frac{\alpha_s(M_Z) - 0.1184}{0.0007} \pm 1.0_{\rm th}$ — at the
  central inputs (identical to this code's: $M_t = 173.1$,
  $\alpha_s = 0.1184$) the bound is $129.4 \pm 1.0$ GeV, i.e. the code
  sits $0.35\sigma_{\rm th}$ below it (enforced by
  `tests/physics/test_boundary.cpp`). This is a consistency check within
  the published theory band, **not** a precision validation: the code
  omits the 2-loop effective potential and 3-loop QCD threshold pieces of
  Degrassi et al. and linearizes the matching. (An earlier version of
  this section anchored on an unverified "$\simeq 128.6$ GeV, Buttazzo
  Fig. 1" — retracted in the 2026-10 referee pass.)
* **Deep negative-coupling regime.** At $(115, 180)$ GeV (classified
  *unstable* by the threshold criterion) the trial action agrees with the
  strict conformal estimate to a few $10^{-4}$ relative, and near the
  boundary at $(134.75, 176.5)$ the two differ by
  $+33\%$ — but the c6 sensitivity test
  (`tests/physics/test_method_flags.cpp`) shows this difference is
  produced by the assumed $c_6 = 1$ Planck operator: with $c_6 = 0$ the
  two actions agree to $-0.9\%$ there and to $-0.08\%$ at the SM point.
  RG improvement of the potential alone moves the trial action by
  $\lesssim 1\%$ everywhere tested.

## 4. What is *not* validated

Stated explicitly to avoid overclaiming:

* **The ansatz itself.** No test validates that the conformal profile is a
  good approximation to the true bounce; it cannot, within this code —
  there is no independent bounce solver to compare against. The trial
  action has **no established relation (bound or otherwise) to the exact
  bounce action**: the bounce is a saddle of the action, so trial actions
  may lie above or below it (see `docs/limitations.md` item 1).
* **The $c_6 = 1$ Planck regularization.** No literature value is
  established; only the existence and sign of the effect are documented.
* **The 4-loop $g_3$ coefficient 2472.28.** Provenance unknown; impact
  measured at $\sim0.1\%$ on benchmark actions by toggling the term.
* **Literature bounce actions.** The code's $S$ values are not expected to
  match full bounce computations and are not marketed as such; they are
  trial-profile estimates on a constrained family with unquantified
  trial-family dependence.

## 6. Resolution convergence of the phase-diagram datasets

The classification datasets exist at five grid spacings (full plane 1.0
and 0.5 GeV; zoom 0.25 and 0.1 GeV; boundary strip 0.05 GeV over
`Mh ∈ [110, 140] × Mt ∈ [162, 182] GeV`), all regenerated after the
2026-10 referee corrections. Every grid point was evaluated by the
pipeline (no interpolation between resolutions); the comparison is
automated in `scripts/check_resolution_convergence.py`. Results at
`Mt = 173 GeV`:

| dataset | grid-level last non-stable Mh | lambda_min = 0 crossing |
|---|---|---|
| full 1.0 GeV | 128.00 | 128.833 |
| full 0.5 GeV | 128.50 | 128.835 |
| zoom 0.25 GeV | 128.75 | 128.836 |
| zoom 0.1 GeV | 128.80 | 128.836 |
| boundary 0.05 GeV | 128.80 | 128.836 |

* The grid-level boundary (last non-stable cell) approaches the true
  crossing from below at the expected one-cell-per-refinement rate; the
  interpolated `lambda_min = 0` level set is resolution-independent to
  **3 mGeV** across all five grids.
* `S_trial` at grid points common to all datasets is bit-identical
  (deterministic pipeline; refinement adds points, never changes values).
* Classification **area fractions** shift by <= 0.25% per halving
  (e.g. full-plane stable 0.3881 -> 0.3897, metastable 0.0571 -> 0.0571,
  unstable 0.1889 -> 0.1892). These fractions are grid-convergence
  diagnostics of the sampled field, not physical predictions over the
  whole plane (the linearized matching limits the physical
  interpretation).
* The ansatz-breakdown point count scales with cell area as expected for
  a fixed-width strip: 347 -> 1,516 points in the full plane (1.0 -> 0.5
  GeV, ratio 4.4 vs 4 expected from h^-2), confirming a physical strip
  width of ~1 grid cell (<= 0.1 GeV).
* The conformal-vs-trial-profile disagreement fraction is stable (1,669
  of 251,001 points, 0.66%, on the 0.5 GeV plane).
* The experimental point is classified **metastable** at every
  resolution.

The phase boundaries and the scientific conclusions are therefore stable
under grid refinement at the stated resolutions.

## 5. Regression tests (`tests/regression/`)

Frozen-value tests guard the pipeline against unintended changes:
RG trajectory and potential tables (`reference/v1.1/*.csv`, regenerated by
`apps/write_reference_tables.cpp` after any intentional physics change) and
single-point benchmark actions (`benchmark_results*.json`). These verify
*reproducibility*, not physics correctness — the physics claims are covered
by §1–§3 above. The current frozen values:
conformal action $= 2120.340020693041$, trial action $= 2168.6895796842282$
(at $c_6 = 1$) at $(125.1, 173.1)$; status 2 (metastable), method flag OK
for both. Frozen after the 2026-10 referee corrections.
