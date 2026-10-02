# Known Limitations

Consolidated, honest list of what this code does *not* do, where it is not
trustworthy, and which results carry which systematics. A reader deciding
how much to trust any number should start here. Items marked
*(affects classification)* can move phase-boundary points; items marked
*(affects rates only)* change absolute action values but not the
stability boundary.

## Method-level approximations

1. **No bounce-equation solve; conformal trial family only.** *(affects
   rates only; trial-family dependence unquantified)* The bounce equation
   $\phi'' + \frac{3}{r}\phi' = \partial V/\partial\phi$ is never
   integrated. The action is evaluated on the Fubini–Lipatov profile
   $\phi_R(r) = \sqrt{2/|\lambda_R|}\;2R/(R^2+r^2)$ with only the scale $R$
   optimized. The result is a **trial-profile action estimate**. **No
   rigorous inequality relating $S_{\rm trial}$ to the exact bounce action
   is established here** — the O(4) bounce is a saddle point of the action
   (one negative mode), not a minimum over general field configurations,
   so a trial action may lie above or below the true bounce action (and
   along the family it can become negative where the ansatz breaks down).
   Coleman, Glaser & Martin's minimum applies only among *solutions of the
   equations of motion*, a class the fixed profile does not belong to.
   Consequently no claim of under-/over-estimated decay rates or
   under-/over-estimated unstable regions follows from this method. No
   independent bounce solver was run, so the size of the trial-family
   dependence is unknown. Only the qualitative phase-diagram topology and
   the stability boundary should be taken over from this code.

2. **High-field potential; false vacuum at the origin.** *(affects
   interpretation)* The tree-level mass term $-\tfrac12 m^2\phi^2$ is
   dropped: $V = \lambda_{\rm eff}(\phi)\phi^4/4$ has its false vacuum at
   $\phi = 0$, not at the electroweak vacuum ($v = 246.22$ GeV). The
   action therefore describes tunneling *from the origin*, used as a proxy
   for the stability of the electroweak vacuum; the true bounce would
   interpolate from $\phi = v$. The barrier shape at $\phi \sim v$ is
   not represented (the barrier maximum is displaced and its height
   inflated relative to the full potential).

3. **Planck-scale operator with assumed $c_6 = 1$.** *(affects rates near
   the boundary, and the classification there)* The potential is augmented
   by $\tfrac{c_6}{6}(\phi/M_{\rm Pl})^2\phi^4$ with the coefficient
   **assumed** to be 1 (an EFT regulator in the spirit of
   Branchina–Messina-type Planck-suppressed operators; sign, size, and
   existence are UV-dependent and no literature provenance exists for the
   specific value). **It dominates the near-boundary action**: with
   $c_6 = 0$ the trial action equals the strict conformal estimate to
   $\lesssim 1\%$ everywhere (e.g. $-0.08\%$ at the SM point, $-0.9\%$ at
   $(134.75, 176.5)$), while $c_6 = 1$ raises the action by up to $+33\%$
   near the boundary. With **negative** $c_6$ the trial action is
   unbounded below along the family ($S_{\rm trial} \simeq -3\times10^{12}$
   at $c_6 = -1$, SM point) — such points are flagged
   `ANSATZ_FAILED`. Earlier drafts attributed the ≈30% near-boundary
   effect to "RG improvement"; that was incorrect. $c_6$ is a free
   parameter of the code (`--c6`).

4. **Prefactor scale in the metastability criterion.** *(affects
   classification near the unstable/metastable boundary)* The criterion
   $S > 4\ln(\mu_1 t_U)$ uses $\mu_1$, the scale where $\lambda_{\rm eff}$
   crosses zero. This is a **conventional choice**, not a determination of
   the fluctuation prefactor: the dimensional scale of the instanton is
   $1/R_\ast$, which at the SM point is $\mu_1/R_\ast \simeq 10^{16}$ GeV.
   Using it instead would **raise** the threshold by
   $4\ln(1/R_\ast) \simeq +47$ at the SM point — about 10% of
   $S_{\rm th} \simeq 484$ — and shift the red/yellow boundary by a few
   tenths of a GeV in $M_t$ (where $\mathrm{d}S/\mathrm{d}M_t \sim
  100$–$200$/GeV). This is the dominant omission relative to
   complete-lifetime treatments (Andreassen, Frost & Schwartz 2018).

5. **Gauge dependence not assessed.** $\lambda_{\rm eff}(\mu)$ and the
   zero-crossing scale $\mu_1$ are gauge-dependent quantities in principle
   (the effective potential itself is gauge-dependent beyond its
   extremum). Since $\mu_1$ enters the metastability threshold, the
   red/yellow classification inherits an unquantified gauge dependence.
   The green (absolute-stability) boundary $\lambda_{\min} = 0$ at the
   Planck scale is the standard, effectively gauge-independent statement.
   The gauge-fixing scheme is not tracked anywhere in the code (the
   matching constants are taken from the published $\overline{\rm MS}$
  /Landau-gauge analysis of Buttazzo et al.).

6. **Flat spacetime, no gravity.** The action is computed in flat
   Euclidean space even though profile amplitudes can reach the Planck
   scale (where the $c_6$ operator and the Planck cap on the running take
   over). Gravitational corrections to the bounce are not included.

## Physics-input limitations

7. **Linearized NNLO matching.** The matching is a fit around
   $(M_h, M_t) = (125.15, 173.34)$ GeV; the $M_h$-dependence enters through
   the tree-level $M_h^2/2v^2$ term, whose non-linearity deviates from the
   published linear fit by up to $8\times10^{-4}$ in $\lambda(M_t)$ at the
   scan edges ($\sim 0.4$ GeV in boundary position). Constants are taken
   as printed in the latest arXiv version of Buttazzo et al. (2013)
   (retrieved via ar5iv); v2/v3 of that paper modified and corrected the
   NNLO $g_2, g_Y$ results, and the version provenance of the individual
   printed digits has not been traced further. The $M_W$-dependence of the
   published fit is **not** implemented ($M_W$ implicitly fixed at its
   central value through the fitted constants); the
   $\alpha_s$-dependence is implemented but hard-wired to the central
   value $0.1184$. The bottom Yukawa is set to the tree-level
   $\sqrt{2}\,M_b/v$ with $M_b = 4.0$ GeV — not an $\overline{\rm MS}$
   $y_b(M_t)$; the effect on $\lambda(M_{\rm Pl})$ is expected to be small
   and has not been quantified.

8. **4-loop $g_3$ running is partial.** The $g_3$ beta function includes
   the pure-QCD nf-dependent 4-loop piece,
   $-\beta_3(n_f{=}6)\,g_3^{10}/(16\pi^2)^4$ with
   $\beta_3 = 2472.28$ (van Ritbergen, Vermaseren & Larin 1997; the sign
   was corrected in 2026-10 — an earlier version had $+2472.28$,
   inconsistent with the code's own convention). The remaining SM 4-loop
   terms (gauge–Yukawa insertions) are absent, so the running is **not**
   a complete 4-loop SM result; the included piece shifts the benchmark
   actions by $\sim 0.1\%$.

9. **Beta-function provenance.** The 3-loop coefficients are numerically
   expanded and hard-coded; spot checks against Buttazzo et al. (2013)
   (matching values, $M_h$-slope, one-loop coefficients) pass, but a
   line-by-line verification of every 3-loop coefficient against a symbolic
   source has not been performed.

## Numerical limitations

10. **Fixed-step RK4 ($\Delta t = 0.1$ in $t = \ln\mu^2$) with linear
    table interpolation.** Measured step-halving convergence: production
    settings reproduce the four-times-finer run to $9.6\times10^{-7}$
    relative in the action (`tests/physics/test_convergence.cpp`). Note
    the table interpolation between stored points is linear
    ($\mathcal{O}(\Delta t^2)$), so the RK4 order does not directly
    describe the accuracy of the queried trajectory.

11. **Quadrature at the round-off floor.** The Simpson potential integral
    ($N = 2048$) is converged to $\sim 10^{-9}$ relative; the integrand has
    an integrable edge singularity at $u \to -1$ (bounce center) whose
    contribution is harmless at this $N$ (validated against the analytic
    conformal limit to $\le 2\times10^{-12}$).

12. **Radius search bracket and unimodality.** The golden-section search
    operates on a fixed bracket of $\pm 8$ in $\ln R$ (about seven decades
    in $R$) around the conformal estimate, assuming unimodality. Inside
    the kinetic-shortcut/fence regions the bracketed function has
    artificial discontinuities, so the returned optimal radius there is
    not meaningful (such points are flagged; their classification is
    dominated by the shortcut overestimate, which exceeds the threshold by
    orders of magnitude either way).

13. **Amplitude/potential mismatch of the trial family.** The profile
    amplitude is set by the raw running coupling
    $\lambda_R = \lambda(\mu_1/R)$ while the potential integrand uses
    $\lambda_{\rm eff} = \lambda + \Delta\lambda_{\rm CW}$. In a
    pure-quartic world the two coincide; here they differ by the CW
    correction, which near $\lambda_R \to 0^-$ makes the amplitude diverge
    while the potential stays finite — this is the origin of the
    non-positive-action breakdown strip along the stability boundary
    (flagged `ANSATZ_FAILED`, never counted as a physical verdict). A more
    self-consistent one-parameter choice (amplitude from
    $\lambda_{\rm eff}(\mu_1/R)$) or a two-parameter trial family
    (independent amplitude and scale) would remove the artifact but
    changes the definition of the estimator; it is left for future work
    rather than patched silently.

## Not done

* No $c_6$ boundary-shift study over both signs on a grid (single-point
  scan only; $c_6$ is available via `--c6`).
* No finite-temperature effects; the calculation is zero-temperature only.
* No full bounce solver (shooting/relaxation) to quantify the ansatz
  systematics of item 1.
* No gauge-dependence study (item 5).
* No uncertainty propagation from experimental inputs
  ($M_t \pm 0.29$ GeV, $M_h \pm 0.11$ GeV PDG 2022) into boundary
  positions; the phase diagram shows central values only.
* No independent cross-check of the RGE trajectory against a third-party
  RG package beyond the published central-value checks of
  `tests/unit/test_matching.cpp`.
