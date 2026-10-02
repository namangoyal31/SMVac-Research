// Literature comparison of the stability classification.
//
// Anchor (verified against the published text): Degrassi et al., JHEP 08
// (2012) 098 [arXiv:1205.6497], introduction eq. (2):
//
//   Mh[GeV] > 129.4 + 1.4 (Mt[GeV] - 173.1)/0.7
//                  - 0.5 (alpha_s(M_Z) - 0.1184)/0.0007  +- 1.0_th
//
// At the central inputs (Mt = 173.1 GeV, alpha_s = 0.1184 -- exactly the
// values used by this code) the bound is Mh_crit = 129.4 +- 1.0 GeV.
// This is the closest like-for-like anchor available: the code uses the
// same inputs and the same lambda_min = 0 definition of absolute
// stability, but omits the 2-loop effective potential and the 3-loop QCD
// threshold correction included in Degrassi et al., and linearizes the
// matching, so agreement at the ~0.1-0.5 GeV level is expected but is not
// a precision validation. The test enforces the published 1-sigma theory
// band and prints the distance explicitly.
//
// The physical point (125.1, 173.1) must be metastable, matching the
// literature consensus that the SM vacuum lifetime exceeds the age of the
// universe.
#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <cmath>
#include <cstdio>

using namespace SMVacuumDecay;

static int failures = 0;

// Locate the absolute-stability boundary in Mh at fixed Mt by bisection:
// status 1 (lambda_min >= 0) for Mh above the boundary, status 2/3 below.
static double stability_boundary_Mh(double Mt, double lo, double hi) {
    for (int it = 0; it < 40; ++it) {
        double mid = 0.5 * (lo + hi);
        StabilityResult r = classify_stability(mid, Mt);
        if (r.status == 1) hi = mid;   // stable -> boundary is below
        else              lo = mid;    // metastable/unstable -> boundary is above
    }
    return 0.5 * (lo + hi);
}

int main() {
    std::printf("=== Physical point classification ===\n");
    StabilityResult sm = classify_stability(125.1, 173.1);
    std::printf("(125.1, 173.1): status = %d (2 = metastable, as in the literature)\n", sm.status);
    if (sm.status != 2) { std::printf("FAIL sm-metastable\n"); ++failures; }

    std::printf("=== Absolute stability boundary at Mt = 173.1 GeV ===\n");
    double Mh_crit = stability_boundary_Mh(173.1, 124.0, 136.0);
    // Degrassi et al. (2012), eq. (2) at Mt = 173.1, alpha_s = 0.1184:
    const double degrassi_central = 129.4;
    const double degrassi_sigma_th = 1.0;
    double distance = Mh_crit - degrassi_central;
    std::printf("measured Mh_crit = %.3f GeV   Degrassi eq.(2): 129.4 +- 1.0 (th) GeV\n",
                Mh_crit);
    std::printf("distance = %+.3f GeV = %.2f sigma_th\n",
                distance, distance / degrassi_sigma_th);
    if (std::abs(distance) > degrassi_sigma_th) { std::printf("FAIL boundary-band\n"); ++failures; }
    std::printf("NOTE: same inputs by construction (alpha_s = 0.1184 hard-wired);\n"
                "this code omits 2-loop potential and 3-loop QCD threshold terms,\n"
                "so this is a consistency check within the published theory band,\n"
                "not a precision validation.\n");

    // Direction sanity: below the boundary the vacuum cannot be absolutely
    // stable, above it lambda_min must be non-negative.
    StabilityResult below = classify_stability(Mh_crit - 1.0, 173.1);
    StabilityResult above = classify_stability(Mh_crit + 1.0, 173.1);
    std::printf("Mh_crit - 1: status %d ; Mh_crit + 1: status %d\n", below.status, above.status);
    if (above.status != 1) { std::printf("FAIL above-boundary-stable\n"); ++failures; }
    if (below.status == 1) { std::printf("FAIL below-boundary-unstable\n"); ++failures; }

    if (failures == 0) { std::printf("PASS\n"); return 0; }
    std::printf("FAIL (%d)\n", failures);
    return 1;
}
