# Known Limitations

Consolidated, honest list of what this code does *not* do, where it is not
trustworthy, and which results carry which systematics. A reader deciding
how much to trust any number should start here. Items marked
*(affects classification)* can move phase-boundary points; items marked
*(affects rates only)* change absolute action values but not the
stability boundary.

## Method-level approximations

1. **No bounce-equation solve; conformal ansatz only.** *(affects rates
   only)* The bounce equation
   $\phi'' + \frac{3}{r}\phi' = \partial V/\partial\phi$ is never
   integrated. The action is evaluated on the Fubini–Lipatov profile
   $\phi_R(r) = \sqrt{2/|\lambda_R|}\;2R/(R^2+r^2)$ with only the scale $R$
   optimized. The result is a variational **upper bound** on the true
   bounce action: decay rates are underestimated, lifetimes overestimated,
   and the computed *unstable* region is correspondingly *underestimated*.
   For SM-like potentials the true bounce is not close to the conformal
   profile (published bounce actions for the physical point are
   substantially smaller than the values reported here); only the
   qualitative topology of the phase diagram and the stability boundary
   should be taken over from this code.

2. **High-field potential; false vacuum at the origin.** *(affects
   interpretation)* The tree-level mass term $-\tfrac12 m^2\phi^2$ is
   dropped: $V = \lambda_{\rm eff}(\phi)\phi^4/4$ has its false vacuum at
   $\phi = 0$, not at the electroweak vacuum ($v = 246.22$ GeV). The
   action therefore describes tunneling *from the origin*, used as a proxy
   for the stability of the electroweak vacuum; the true bounce would
   interpolate from $\phi = v$. The barrier shape at $\phi \sim v$ is
   not represented (the barrier maximum is displaced and its height
   inflated relative to the full potential).

3. **Planck-scale regularization with $c_6 = 1$.** *(affects rates near
   the boundary)* The potential is augmented by
   $\tfrac{c_6}{6}(\phi/M_{\rm Pl})^2\phi^4$ with the coefficient
   **assumed** to be 1. No literature provenance was established for this
   value; it dominates the action for bubbles reaching the Planck scale.
   A sensitivity study over $c_6 \in [0.1, 10]$ is not yet performed and
   is the single most valuable physics extension of this code (earlier
   internal notes estimated boundary shifts of order $0.1$–$1$ GeV from
   varying $c_6$).

4. **Prefactor scale in the metastability criterion.** *(affects
   classification near the unstable/metastable boundary)* The criterion
   $S > 4\ln(\mu_1 t_U)$ uses $\mu_1$, the scale where $\lambda_{\rm eff}$
   crosses zero. The prefactor scale of the decay rate is more properly
   $\sim 1/R_\ast$ (the inverse optimal bubble radius); neglecting the
   radius dependence of the prefactor is the dominant omission relative to
   complete-lifetime treatments (Andreassen, Frost & Schwartz 2018).
   Using $1/R_\ast$ instead of $\mu_1$ would lower the threshold by
   $4\ln(\mu_1 R_\ast)$ — a shift comparable to a few grid spacings of the
   phase boundary.

5. **Flat spacetime, no gravity.** The bounce is computed in flat
   Euclidean space even though profile amplitudes can reach the Planck
   scale (where the $c_6$ operator and the Planck cap on the running take
   over). Gravitational corrections to the bounce are not included.

## Physics-input limitations

6. **Linearized NNLO matching.** The matching is a fit around
   $(M_h, M_t) = (125.15, 173.34)$ GeV; the $M_h$-dependence enters through
   the tree-level $M_h^2/2v^2$ term, whose non-linearity deviates from the
   published linear fit by up to $8\times10^{-4}$ in $\lambda(M_t)$ at the
   scan edges ($\sim 0.4$ GeV in boundary position). The $\alpha_s$- and
   $M_W$-dependences of the matching are implemented but hard-wired to
   central values.

7. **Unverified 4-loop $g_3$ term.** The $g_3$ beta function includes a
   4-loop contribution with numeric coefficient $2472.28$ whose provenance
   could not be established from the code or the cited literature. It is
   retained (removing it would silently change the physics); toggling it
   shifts the benchmark actions by $\sim 0.1\%$.

8. **Beta-function provenance.** The 3-loop coefficients are numerically
   expanded and hard-coded; spot checks against Buttazzo et al. (2013)
   (matching values, $M_h$-slope, one-loop coefficients) pass, but a
   line-by-line verification of every 3-loop coefficient against a symbolic
   source has not been performed.

## Numerical limitations

9. **Fixed-step RK4 ($\Delta t = 0.1$ in $t = \ln\mu^2$).** Measured
   step-halving convergence: the production settings reproduce the
   four-times-finer run to $9.6\times10^{-7}$ relative in the action
   (`tests/physics/test_convergence.cpp`). Not adaptive; the (removed)
   adaptive stepper of the original code was never used by the pipeline.

10. **Quadrature at the round-off floor.** The Simpson potential integral
    ($N = 2048$) is converged to $\sim 10^{-9}$ relative; the integrand has
    an integrable edge singularity at $u \to -1$ (bounce center) whose
    contribution is harmless at this $N$ (validated against the analytic
    conformal limit to $\le 2\times10^{-12}$).

11. **Radius search bracket and unimodality.** The golden-section search
    operates on a fixed 16-decade bracket in $\ln R$ around the conformal
    estimate, assuming unimodality. No pathology of this assumption has
    been observed across the committed scan grids, but it is not proven.

12. **Ansatz breakdown near the stability boundary.** Immediately adjacent
    to the $\lambda_{\min} = 0$ boundary, isolated grid points yield
    **non-positive ansatz actions** (the profile amplitude
    $\sqrt{2/|\lambda_R|}$ diverges as $\lambda_R \to 0^-$ while the
    potential integral uses $\lambda_{\rm eff} \neq \lambda_R$). In the
    committed 0.25 GeV SM-window scan, 62 of 2,324 "unstable" points
    (2.7%) have $S \le 0$ (full plane: 518); they form a one-grid-cell
    strip along the stability boundary and are drawn as a separate
    category in the phase diagram. They are an artifact of the method, not
    physical instabilities. Two further artifact classes exist near the
    boundary and are flagged rather than hidden: **kinetic-shortcut
    points** (61 in the SM window, 98 in the full plane), where the action
    is returned without the potential integral and is only a rough
    overestimate, and **fence points** (11 in the SM window and full
    plane), where no radius in the search bracket gives
    $\lambda_R < 0$, so the "action" is the $10^{100}$ sentinel rather
    than a calculation (the classification as metastable is still correct —
    the true action is large). All three classes are exported in the CSV
    columns and drawn explicitly in figures 02, 04 and 07. The
    kinetic-shortcut approximation also overestimates such actions by up
    to a factor $\sim 2$ — irrelevant for the classification since those
    actions exceed the threshold by orders of magnitude either way.

13. **Interpolation of the RG table.** Couplings are linearly interpolated
    between table points ($\Delta t = 0.1$); the associated error is
    covered by the step-size study but has not been separated from the RK4
    truncation error.

## Not done

* No $c_6$ sensitivity scan (item 3) — highest-value physics extension.
* No finite-temperature effects; the calculation is zero-temperature only.
* No full bounce solver (shooting/relaxation) to quantify the ansatz
  systematics of item 1.
* No uncertainty propagation from experimental inputs
  ($M_t \pm 0.29$ GeV, $M_h \pm 0.11$ GeV PDG 2022) into boundary
  positions; the phase diagram shows the central values only, with the
  PDG ellipses drawn for orientation.
* No comparison against an independent RGE package (e.g. one-loop checks
  against a third-party tool) beyond the published central-value checks of
  `tests/unit/test_matching.cpp`.
