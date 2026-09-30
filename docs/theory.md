# Theory

This document describes the physical problem, the theoretical setup, every
equation implemented in the code, the conventions, and the approximations
that make the calculation tractable. The mapping from these equations to the
implementation is given in [numerical-method.md](numerical-method.md); the
validation status of each ingredient is summarized in
[validation.md](validation.md).

---

## 1. Physical problem

The measured Higgs and top masses place the Standard Model (SM) close to the
boundary between *absolute stability* and *metastability* of the electroweak
vacuum. At zero temperature the decay rate of a metastable vacuum is
controlled by the Euclidean bounce action $S$ (the O(4)-symmetric instanton
of Coleman), and the SM is metastable if its lifetime exceeds the age of the
universe. Whether it does depends on the running of the Higgs quartic
coupling $\lambda(\mu)$, which is driven negative at large field values by
the top Yukawa — but only barely so for the central measured masses.

This repository computes, in the $(M_h, M_t)$ plane:

* the **absolute stability boundary**: where the RG-improved quartic coupling
  stays positive all the way to the Planck scale ($\lambda_{\min} = 0$);
* an estimate of the **decay action** $S$ in the metastable region, computed
  on the Fubini–Lipatov (conformal) profile with RG-improved running;
* a **metastable / unstable classification** from the age-of-the-universe
  criterion applied to that action.

The two-step nature of the method must be kept in mind when interpreting
results: the *stability boundary* follows directly from the RG potential and
is robust; the *action in the metastable region* is computed for a fixed
(conformal) profile family and is a variational upper bound on the true
bounce action. See §6 and [limitations.md](limitations.md).

## 2. Setup and conventions

* **Units.** Natural units $\hbar = c = 1$; all dimensional quantities in
  GeV. Actions are dimensionless.
* **Field.** $\phi$ is the real, radially parametrized SM Higgs background
  ($\phi = \langle h^0\rangle$). The electroweak vacuum sits at
  $v = 246.22\ \mathrm{GeV}$.
* **Potential normalization.** $V \supset \lambda \phi^4/4$, so that at tree
  level $M_h^2 = 2\lambda v^2$.
* **Hypercharge.** $g_1 = \sqrt{5/3}\, g_Y$ (GUT-normalized).
* **RG time.** The evolution variable is $t = \ln \mu^2$, and *every* beta
  function in the code is implemented as $\mathrm{d}X/\mathrm{d}t =
  \mathrm{d}X/\mathrm{d}\ln\mu^2 = \tfrac12\, \mathrm{d}X/\mathrm{d}\ln\mu$.
  The rates are therefore **half** the standard $\mathrm{d}/\mathrm{d}\ln\mu$
  coefficients; the factor is absorbed consistently in the RK4 stepper. This
  convention is self-consistent but unusual — it is documented here because
  it is easy to misread as a factor-2 bug.
* **RG improvement.** Couplings are evaluated at $\mu = \phi$ (the field
  value), which resums the leading $\ln(\phi/v)$ terms of the potential.
* **Planck scale.** $M_{\mathrm{Pl}} = 1.22\times10^{19}$ GeV
  (non-reduced). The RGEs are integrated from $\mu = M_t$ to
  $\mu = M_{\mathrm{Pl}}$; couplings are **frozen** above $M_{\mathrm{Pl}}$.
  Gravity is not included in the action (flat-space bounce).

## 3. Matching conditions at the top scale

The $\overline{\mathrm{MS}}$ couplings at $\mu = M_t$ are set by the NNLO
electroweak matching conditions of Buttazzo et al. (2013), eqs. (55)–(60),
linearized around the central masses $M_t = 173.34$ GeV,
$M_h = 125.15$ GeV. The implementation expresses the $M_h$-dependence
through the tree-level relation plus a constant offset:

$$
\lambda(M_t) \;=\; \frac{M_h^2}{2v^2} \;-\; 0.00313 \;-\; 0.00004\,(M_t - 173.34),
$$

$$
y_t(M_t) \;=\; \frac{\sqrt{2}\,M_t}{v} \;-\; 0.0587 \;-\; 0.00042\,
\frac{\alpha_s(M_Z) - 0.1184}{0.0007},
$$

$$
g_1(M_t) = 0.462458 + 0.000142\,(M_t - 173.34), \qquad
g_2(M_t) = 0.64779 + 0.00004\,(M_t - 173.34),
$$

$$
g_3(M_t) = 1.1666 + 0.00314\,\frac{\alpha_s(M_Z) - 0.1184}{0.0007}
          - 0.00046\,(M_t - 173.34).
$$

Notes and caveats (see also `docs/audit-notes.md`, finding 16):

* The tree-level $M_h^2/(2v^2)$ term reproduces the fitted NNLO slope
  $\mathrm{d}\lambda/\mathrm{d}M_h = 0.00206\ \mathrm{GeV}^{-1}$; its
  non-linearity deviates from the linear fit by up to $8\times10^{-4}$ in
  $\lambda(M_t)$ at the edges of the scan range ($\sim 0.4$ GeV in boundary
  position).
* At the central masses the implementation gives
  $\lambda(M_t) = 0.126025$ vs. the published $0.12604$ (difference
  $1.5\times10^{-5}$), and $y_t(M_t) = 0.936888$ vs. $0.93690$. These checks
  are automated in `tests/unit/test_matching.cpp`.
* $\alpha_s(M_Z)$ is fixed to its central value $0.1184$; the
  $\alpha_s$-dependence of the matching is implemented but never activated.
* The bottom and tau Yukawas are set at $\mu = M_t$ to
  $y_b = \sqrt{2}\,M_b/v$ with $M_b = 4$ GeV and $y_\tau = \sqrt{2}\,M_\tau/v$
  with $M_\tau = 1.777$ GeV; their running is included only in the beta
  functions (no threshold matching).

## 4. Running and the effective potential

### 4.1 RGEs

The 3-loop SM beta functions for $(g_1, g_2, g_3,\, y_t, y_b, y_\tau,\,
\lambda)$ are implemented exactly as compiled in Buttazzo et al. (2013)
Appendix (numerically expanded, in the squared-coupling variable $X = x^2$
with the $\mathrm{d}/\mathrm{d}\ln\mu^2$ convention of §2), plus one
additional 4-loop pure-gauge $g_3$ term whose coefficient ($2472.28$) could
not be traced to a published source; its measured impact on the benchmark
actions is $\sim 0.1\%$ (`docs/audit-notes.md`, finding 13). The system is

$$
\frac{\mathrm{d}X_i}{\mathrm{d}t} = \beta_i^{(1)} + \frac{\beta_i^{(2)}}{(16\pi^2)} +
\frac{\beta_i^{(3)}}{(16\pi^2)^2} + \dots ,
\qquad t = \ln\mu^2 ,
$$

with the loop factors normalized as in the code. The trajectory is
integrated by fixed-step classical RK4 with $\Delta t = 0.1$ from
$t_0 = \ln M_t^2$ up to $t_{\mathrm{Pl}} = \ln M_{\mathrm{Pl}}^2$, storing a
table that is later queried (with linear interpolation in $t$) at arbitrary
scales. Step-size convergence is quantified in
[validation.md](validation.md).

### 4.2 RG-improved effective potential

For $\phi \gg v$ the potential is taken in the high-field form

$$
V(\phi) \;=\; \frac{1}{4}\,\lambda_{\mathrm{eff}}(\phi)\,\phi^4 ,
\qquad
\lambda_{\mathrm{eff}}(\phi) \;=\; \lambda(\mu{=}\phi) + \Delta\lambda_{\mathrm{CW}},
$$

with the 1-loop Coleman–Weinberg correction from the top quark and the
$W/Z$ gauge bosons, written with field-dependent masses
$m_t^2 = \tfrac12 y_t^2\phi^2$, $m_W^2 = \tfrac14 g_2^2\phi^2$,
$m_Z^2 = \tfrac14 (g_1^2+g_2^2)\phi^2$:

$$
\Delta\lambda_{\mathrm{CW}} = \frac{1}{16\pi^2}\left[
-3 y_t^4\left(\ln\frac{y_t^2}{2} - \frac32\right)
+ \frac{3}{8} g_2^4\left(\ln\frac{g_2^2}{4} - \frac56\right)
+ \frac{3}{16}(g_1^2{+}g_2^2)^2\left(\ln\frac{g_1^2+g_2^2}{4} - \frac56\right)
\right].
$$

(These are the $n_i m_i^4(\ln(m_i^2/\mu^2) - c_i)/(64\pi^2)$ terms with
$\mu = \phi$, $n_t = -12$, $c_t = 3/2$, $n_W = 6$, $n_Z = 3$, $c_V = 5/6$,
divided by $\phi^4/4$; the Goldstone and Higgs loops carry no $\phi^4$
weight and are omitted — negligible for $\phi \gg v$.)

**Deliberate approximation.** The tree-level mass term $-\tfrac12 m^2\phi^2$
is dropped. In this high-field potential the origin $\phi = 0$ is the false
vacuum; the true electroweak vacuum structure (minimum at $v$, barrier at
$\phi\sim 2v$) is absent. This is the price of the conformal ansatz and is
harmless for the deep-field part of the bounce (where $|\lambda_{\rm eff}|\phi^4/4$
dominates the action), but it means the computed action is *not* the action
of a bounce interpolating between the true EW vacuum and the deep-field
region. The stability *boundary* is unaffected: it depends only on whether
$\lambda_{\rm eff}$ runs negative.

## 5. Bounce action with the conformal profile

### 5.1 Bounce formalism

At zero temperature the decay rate per unit volume of a metastable vacuum is

$$
\Gamma/V \;\simeq\; \mu^4\, e^{-S},
$$

where $S$ is the Euclidean action of the O(4)-symmetric bounce,

$$
S[\phi] \;=\; 2\pi^2 \int_0^\infty \mathrm{d}r\; r^3
\left[\frac12 \left(\frac{\mathrm{d}\phi}{\mathrm{d}r}\right)^2 + V(\phi)\right],
$$

and the bounce profile solves

$$
\frac{\mathrm{d}^2\phi}{\mathrm{d}r^2} + \frac{3}{r}\frac{\mathrm{d}\phi}{\mathrm{d}r}
= \frac{\partial V}{\partial \phi},
\qquad
\phi'(0) = 0, \quad \phi(\infty) = 0 .
$$

The true bounce minimizes $S$ over *all* profiles. **This code never solves
the bounce equation.** Instead it evaluates $S$ on the Fubini–Lipatov
(conformal) profile family and minimizes over its single scale parameter.
The result is therefore a **variational upper bound**: $S_{\rm ansatz} \ge
S_{\rm true}$, i.e. decay rates are underestimated and lifetimes
overestimated. For the pure-quartic potential the family contains the exact
solution and the bound is saturated.

### 5.2 The Fubini–Lipatov profile

For the exactly scale-invariant potential
$V = -|\lambda|\phi^4/4$ the bounce equation has the analytic solution

$$
\phi_R(r) \;=\; \sqrt{\frac{2}{|\lambda|}}\;\frac{2R}{R^2 + r^2},
\qquad \phi_R(0) = \frac{2\sqrt{2/|\lambda|}}{R},
$$

with $R$ an arbitrary scale (a consequence of classical scale invariance).
For the RG-improved potential the code uses this same profile with the
amplitude set by the running coupling at the inverse bubble scale,

$$
\lambda_R \;=\; \lambda(\mu = \mu_{\rm inst}/R), \qquad
A = \sqrt{2/|\lambda_R|}.
$$

The kinetic term of the fixed profile is analytic for *any* potential:

$$
S_{\rm kin} \;=\; 2\pi^2\!\int_0^\infty\! \mathrm{d}r\, r^3\,\frac{\phi_R'^2}{2}
= \frac{8\pi^2}{3} A^2 \;=\; \frac{16\pi^2}{3\,|\lambda_R|}.
$$

The potential term is evaluated numerically:

$$
S_{\rm pot} \;=\; 2\pi^2\!\int_0^\infty\! \mathrm{d}r\, r^3\, V(\phi_R(r))
\;=\; 2\pi^2 R^4 \int_{-1}^{1}\mathrm{d}u\;\,
\frac{\alpha}{1-u^2}\,\left(\frac{1+u}{1-u}\right)^{2\alpha}
V\!\left(\phi_R(r(u))\right),
$$

under the compactifying substitution $r = R\,e^{\alpha x}$, $u = \tanh x$
($\alpha = 4$ in the code). The total ansatz action is
$S(R) = S_{\rm kin} + S_{\rm pot}$, and the code minimizes $S(R)$ over
$\ln R$ by golden-section search in an 8-decade bracket around the naive
conformal estimate $R \sim \mu_{\rm inst}/\mu(\lambda_{\min})$.

### 5.3 Planck-scale regularization ($c_6$ operator)

As $\phi \to M_{\mathrm{Pl}}$ the running of $\lambda$ is extrapolated beyond
its domain of validity and the pure-quartic ansatz becomes ill-behaved. The
code therefore augments the potential with a Planck-suppressed dimension-6
operator,

$$
V(\phi) \;\supset\; \frac{c_6}{6}\left(\frac{\phi}{M_{\mathrm{Pl}}}\right)^2 \phi^4 ,
\qquad c_6 = 1 ,
$$

which dominates the quartic term as $\phi \to M_{\mathrm{Pl}}$ and
regularizes the action of large bubbles. **The specific value $c_6 = 1$ is
an assumption of the original implementation**: no literature provenance was
established for it, it is hard-coded as `c6_planck_suppressed` in
`src/CanonicalBounce.cpp`, and the sensitivity of results to $c_6$ is
listed under known limitations. Note that $c_6 > 0$ *raises* the potential
at large field values and hence *raises* $S$; it can only make the vacuum
look more stable.

### 5.4 Strict conformal estimate

For orientation the code also reports the strict Fubini–Lipatov estimate
formed from the minimum of the running coupling along the trajectory,

$$
S_{\rm approx} \;=\; \frac{8\pi^2}{3\,|\lambda_{\min}|},
\qquad \lambda_{\min} = \min_{\mu}\lambda_{\rm eff}(\mu),
$$

which is the exact bounce action in the pure-quartic limit and agrees with
the optimized ansatz action deep in the metastable region (validated to
$\lesssim 0.1\%$, see [validation.md](validation.md)). Near the stability
boundary the two differ by tens of percent — that difference is the
motivation for the RG-improved evaluation.

## 6. Metastability criterion

With the ansatz action $S_\ast = \min_R S(R)$ and the dimensional prefactor
scale $\mu$, the probability for a decay inside the past light cone of a
universe of age $t_U$ scales as $(\mu\, t_U)^4 e^{-S}$. Requiring it to stay
below unity gives the **metastability threshold**

$$
S \;>\; S_{\rm th} \;=\; 4 \ln\!\left(\mu\, t_U\right).
$$

The code evaluates this with:

* $\mu = \mu_1$, the scale at which $\lambda_{\rm eff}$ first crosses zero
  (the only physical scale in the problem besides the bubble radius; the
  choice of prefactor scale is a systematic uncertainty — using
  $\mu \sim 1/R_\ast$ instead lowers $S_{\rm th}$ by $4\ln(\mu_1 R_\ast)$,
  shifting the unstable/metastable boundary; see limitations);
* $t_U = 10$ Gyr (a conservative lower bound; the current 13.8 Gyr estimate
  would raise $S_{\rm th}$ by only $\simeq 1.3$), implemented via named
  constants in `include/SMVacuumDecay/RGE.hpp`.

Classification outcome:

| status | meaning |
|--------|---------|
| 1 | stable: $\lambda_{\rm eff}(\mu)$ never crosses zero up to $M_{\rm Pl}$ |
| 2 | metastable: $\lambda_{\min} < 0$ and $S_\ast > S_{\rm th}$ |
| 3 | unstable: $\lambda_{\min} < 0$ and $S_\ast \le S_{\rm th}$ |
| 4 | perturbativity lost before $M_{\rm Pl}$ ($|\lambda| > 4\pi$, $y_t > 4\pi$, or non-finite couplings) |

Points of type 3 include regions where the ansatz action becomes **negative**
(the potential integral overwhelms the kinetic term). Negative actions are a
documented breakdown of the conformal ansatz far outside its domain, not a
numerical error; they always classify as unstable, and the exported
$S_{\rm kin}/|S_{\rm pot}|$ ratio allows them to be identified.

## 7. Interpretation of outputs

* `S_exact` — the optimized ansatz action $S_\ast$ (the name is historical;
  it is *not* the exact bounce action). Dimensionless.
* `S_approx` — the strict conformal estimate $8\pi^2/(3|\lambda_{\min}|)$.
* `S_kinetic`, `S_potential` — the closed-form kinetic term and the
  Simpson-integrated potential term at the optimal radius; their ratio
  $S_{\rm kin}/|S_{\rm pot}| = 2$ in the pure-quartic limit.
* `S_threshold` — the metastability threshold $4\ln(\mu_1 t_U)$.
* `mu_inst` — the $\lambda_{\rm eff}$ zero-crossing scale $\mu_1$ (GeV).
* `lambda_min` — the minimum of $\lambda_{\rm eff}(\mu)$ over the running
  range.

## 8. Known limitations

See [limitations.md](limitations.md) for the consolidated list — most
importantly: the conformal ansatz (no bounce-equation solve), the
high-field potential (false vacuum at the origin), the assumed
$c_6 = 1$, the prefactor-scale choice in the threshold, the linearized
matching, and the unverified 4-loop $g_3$ term.

## References

See [references.md](references.md) for the full, verified bibliography.
