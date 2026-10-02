# Results

What the corrected calculation actually establishes, with the boundaries
between established, estimated, and unknown made explicit. All numbers are
reproducible via the commands in the README and `data/README.md`.

## 1. The absolute-stability boundary (the robust result)

The RG-improved quartic coupling `λ_eff(μ)` (1-loop Coleman–Weinberg
over 3-loop running, NNLO matching at `μ = Mt`) stays non-negative up to
the Planck scale above, at `αs = 0.1184`:

| Mt [GeV] | `Mh_crit` [GeV] (this code) |
|---:|---:|
| 173.1  | 129.05 |
| 173.34 | 129.55 |

Consistency anchor: Degrassi et al. (2012), eq. (2) —
`Mh > 129.4 + 1.4(Mt−173.1)/0.7 − 0.5(αs−0.1184)/0.0007 ± 1.0_th` — at the
same inputs gives `129.4 ± 1.0 GeV`; the code sits 0.35σ_th below it. This
is a consistency check, not a precision validation: the code omits the
2-loop effective potential and 3-loop QCD threshold pieces of that
analysis and linearizes the matching. The boundary depends only on
`λ_min = 0` and is independent of the trial-profile ansatz.

## 2. The trial-profile action (an estimate, not a bounce action)

At the benchmark point (125.1, 173.1) GeV the vacuum is metastable:
`λ_eff` crosses zero at `μ1 = 7.6e10 GeV` and the optimized trial-profile
action is `S_trial = 2168.7` (c6 = 1), far above the threshold
`S_th = 484` from the Coleman criterion (`Γ/V ~ t_U^-4`, t_U = 10 Gyr).
**`S_trial` is the Euclidean action of a restricted conformal trial
family, minimized over the profile scale; no inequality relating it to the
exact bounce action is established** (the bounce is a saddle point; see
`docs/theory.md` §5.1). No independent bounce computation was performed,
so the accuracy of the estimate is unknown.

## 3. Inside the trial family, RG improvement ≈ negligible; the assumed Planck operator dominates

With the Planck-suppressed operator switched off (`c6 = 0`), the
trial-profile action and the strict conformal estimate
`S_conformal = 8π²/(3|λ_min|)` agree to ≲1% at every point tested
(−0.08% at the SM point; −0.9% at (134.75, 176.5) GeV near the boundary).
The assumed `c6 = 1` operator raises the action by up to +33% near the
boundary. Scan at the SM point (S_trial):

| c6 | −1 | −0.1 | 0 | +0.1 | +1 | +10 |
|---|---|---|---|---|---|---|
| S_trial | −2.8e12* | −2.8e11* | 2118.7 | 2142.4 | 2168.7 | 2208.9 |

\* negative `c6` makes the trial action unbounded below along the family
(flagged `ANSATZ_FAILED`). The c6 coefficient is an assumption; nothing in
this repository predicts it. Earlier drafts quoted the near-boundary
enhancement as an "≈30% RG-improvement effect" — that was incorrect, as
the c6 = 0 row shows.

## 4. Main limitations

No bounce-equation solve (trial-family dependence unquantified); no
bound relation to the true bounce action; high-field potential with the
false vacuum at the origin; assumed c6; prefactor-scale convention
(μ1 instead of 1/R*, worth ≈ +47 in S_th at the SM point); gauge
dependence of μ1 unassessed; linearized NNLO matching with αs, M_W, Mb
fixed (αs-dependence implemented but inert); incomplete 4-loop g3
running; ansatz-breakdown/fence/shortcut points flagged per point and
excluded from physical classifications. Full list:
`docs/limitations.md`.
