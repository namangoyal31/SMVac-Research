# References

Only sources that were actually consulted for, or are directly implemented
or checked against, this codebase are listed. No citation is decorative.

## Method: vacuum decay and the conformal instanton

1. **S. Coleman**, *The Uses of the Instanton*, Proc. 1977 Int. School of
   Subnuclear Physics "Ettore Majorana" (Erice), reprinted in *Aspects of
   Symmetry* (Cambridge University Press, 1985) 265–350.
   — Bounce formalism, decay rate $\Gamma/V$, the thin-wall limit, and the
   dimensional prefactor argument used for the metastability threshold.

2. **S. Fubini**, *A New Approach to Conformal Invariant Field Theories*,
   Nuovo Cimento A **34** (1976) 521.
   — The conformal instanton $\phi(r) \propto 1/(1 + r^2/R^2)$ for a pure
   quartic potential; the exact action $S = 8\pi^2/(3|\lambda|)$.

3. **L. N. Lipatov**, *Divergence of the Perturbation Theory Series and
   Pseudoparticles*, Sov. Phys. JETP **45** (1977) 216.
   — Semiclassical evaluation of the large-order behavior of perturbation
   theory for anharmonic potentials; the same conformal saddle.

4. **C. G. Callan, Jr. and S. Coleman**, *The Fate of the False Vacuum. 2.
   First Quantum Corrections*, Phys. Rev. D **16** (1977) 1762.
   — Fluctuation prefactor and the interpretation of the decay rate used in
   the threshold criterion.

## Standard Model vacuum stability (the physics context and benchmarks)

5. **G. Isidori, G. Ridolfi and A. Strumia**, *On the Metastability of the
   Standard Model Vacuum*, Nucl. Phys. B **609** (2001) 387
   [arXiv:hep-ph/0104016].
   — Canonical RG-improved metastability analysis; the identification of
   $\lambda_{\min}$ as the control parameter for the decay action.

6. **G. Degrassi, S. Di Vita, J. Elías-Miró, J. R. Espinosa, G. F. Giudice,
   G. Isidori and A. Strumia**, *Higgs Mass and Vacuum Stability in the
   Standard Model at NNLO*, JHEP **08** (2012) 098
   [arXiv:1205.6497].
   — NNLO stability analysis; reference for the structure of the matching
   conditions.

7. **D. Buttazzo, G. Degrassi, P. P. Giardino, G. F. Giudice, F. Sala,
   A. Salvio and A. Strumia**, *Investigating the Near-criticality of the
   Higgs Boson*, JHEP **12** (2013) 089 [arXiv:1307.3536].
   — **Primary reference for this codebase**: the NNLO matching conditions
   (eqs. 55–60) implemented in `get_nnlo_matching`, the 3-loop SM beta
   functions implemented in `src/RGE.cpp`, and the stability/metastability
   phase structure in the $(M_h, M_t)$ plane that the phase diagram of this
   repository reproduces qualitatively.

8. **A. Andreassen, W. Frost and M. D. Schwartz**, *Scale Invariant
   Instantons and the Complete Lifetime of the Standard Model*, Phys. Rev. D
   **97** (2018) 056006 [arXiv:1707.08123].
   — State-of-the-art bounce computation including the prefactor; the basis
   for the caveat that neglecting the prefactor's scale dependence (as this
   code does) shifts lifetime estimates, and for the observation that the
   conformal profile is the correct saddle in the pure-quartic limit.

## Numerical methods

9. **W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery**,
   *Numerical Recipes*, 3rd ed. (Cambridge University Press, 2007),
   chs. 16–17.
   — Classical RK4 fixed-step integration, composite Simpson quadrature and
   golden-section search (the three generic algorithms used here, in
   `include/SMVacuumDecay/numerics/`).

## Constants

10. **Particle Data Group** (R. L. Workman *et al.*), *Review of Particle
    Physics*, Prog. Theor. Exp. Phys. **2022** (2022) 083C01, and later
    editions.
    — $\alpha_s(M_Z) = 0.1184$, $M_\tau = 1.777$ GeV, $v = 246.22$ GeV,
    $\hbar$ (seconds↔GeV$^{-1}$ conversion for the universe-age criterion),
    $M_{\rm Pl} = 1.22\times10^{19}$ GeV (non-reduced, gravitational
   units convention).

---

### Provenance note

The beta-function coefficients and matching constants enter the code as
numerically expanded, hard-coded expressions whose line-by-line provenance
could not be re-established for every coefficient during the audit — most
notably the 4-loop $g_3$ term with coefficient `2472.28`
(`docs/audit-notes.md`, finding 13). Where a check was possible (matching
central values, the $M_h$-slope of $\lambda(M_t)$, one-loop beta-function
coefficients) the code agrees with the sources above as documented in
`docs/validation.md`; where it was not possible, the term is flagged rather
than silently trusted or removed.
