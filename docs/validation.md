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
  $S_\ast = 2103.1111,\ 2103.1046,\ 2103.1030,\ 2103.1026$: successive
  differences shrink by a factor $\simeq 4$ per halving, and the production
  step $\Delta t = 0.1$ agrees with the four-times-finer run to
  $9.6\times10^{-7}$ relative. The RG evolution is therefore not a
  significant source of error at production settings.
* **Quadrature.** At the optimal radius, $S(R_\ast)$ for $N \in \{256,
  1024, 4096, 16384\}$ Simpson points agrees to $4.3\times10^{-9}$ relative
  between $N = 1024$ and $N = 16384$. The differences between successive $N$
  are non-monotone at the $\sim10^{-9}$ relative level: the production
  setting $N = 2048$ already sits at the round-off floor (the potential
  integral involves large cancellations), so no meaningful order estimate
  is possible there. The theoretical $\mathcal{O}(N^{-4})$ Simpson
  convergence is verified separately on smooth integrands in
  `tests/unit/test_numerics.cpp` (observed order 4.00).
* **Radius minimization.** The golden-section search agrees with a 50,001-point
  grid scan of the same 16-decade bracket to $1.2\times10^{-10}$ relative.

## 3. Literature comparison (`tests/physics/test_boundary.cpp`)

* **Physical point.** $(M_h, M_t) = (125.1, 173.1)$ GeV is classified
  **metastable** (status 2), matching the consensus that the SM vacuum
  lifetime exceeds the age of the universe.
* **Absolute stability boundary.** Bisection at $M_t = 173.1$ GeV locates
  $\lambda_{\min} = 0$ at $M_h^{\rm crit} = 129.17$ GeV, within $0.6$ GeV of
  the published NNLO value $M_h^{\rm crit} \simeq 128.6$ GeV at
  $M_t \simeq 173.3$ GeV (Buttazzo et al. 2013, Fig. 1). The residual
  $\sim0.6$ GeV shift is consistent with the linearized matching and the
  high-field potential used here (the test enforces a $\pm2.5$ GeV band).
  For reference, the pre-audit matching (with the untraceable $g_2(M_t)$)
  placed the boundary at $\simeq 128.8$ GeV but overestimated the actions
  away from the boundary by 2–24%.
* **Deep metastability limit.** At $(115, 180)$ GeV the optimized ansatz
  action agrees with the strict conformal estimate to $2.1\times10^{-4}$
  relative (running corrections are negligible deep in the metastable
  region), and near the boundary at $(134.75, 176.5)$ the two differ by
  $+29.7\%$ — quantifying when the running-coupling evaluation matters.

## 4. What is *not* validated

Stated explicitly to avoid overclaiming:

* **The ansatz itself.** No test validates that the conformal profile is a
  good approximation to the true bounce; it cannot, within this code —
  there is no independent bounce solver to compare against. The action is a
  variational upper bound, and absolute lifetime predictions inherit that
  systematic (see `docs/limitations.md`).
* **The $c_6 = 1$ Planck regularization.** No literature value is
  established; only the existence and sign of the effect are documented.
* **The 4-loop $g_3$ coefficient 2472.28.** Provenance unknown; impact
  measured at $\sim0.1\%$ on benchmark actions by toggling the term.
* **Literature bounce actions.** The code's $S$ values are not expected to
  match full bounce computations (e.g. $\sim 200$–$400$ in the literature
  for SM-like points at the true saddle) and are not marketed as such; they
  are upper bounds computed on the constrained profile family.

## 5. Regression tests (`tests/regression/`)

Frozen-value tests guard the pipeline against unintended changes:
RG trajectory and potential tables (`reference/v1.1/*.csv`, regenerated by
`apps/write_reference_tables.cpp` after any intentional physics change) and
single-point benchmark actions (`benchmark_results*.json`). These verify
*reproducibility*, not physics correctness — the physics claims are covered
by §1–§3 above. The current frozen values:
FL action $= 2051.1373669116429$, ansatz action $= 2103.1046353416655$ at
$(125.1, 173.1)$; status 2 for both.
