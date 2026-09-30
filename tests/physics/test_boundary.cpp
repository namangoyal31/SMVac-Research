// Literature comparison of the stability classification.
//
// The absolute-stability boundary of the SM (where lambda_eff stays positive
// up to the Planck scale) is known from dedicated NNLO studies (Buttazzo et
// al. 2013, Fig. 1) to lie near Mh_crit ~ 128.6 GeV at Mt ~ 173.3 GeV, with
// a theoretical spread of a few tenths of a GeV between state-of-the-art
// implementations. This code uses a 3-loop RGE with a linearized matching
// and the high-field potential, so only a band comparison is appropriate:
// the measured boundary must fall within +/- 2.5 GeV of the published
// central value. The physical point (125.1, 173.1) must be metastable,
// matching the literature consensus that the SM vacuum lifetime exceeds the
// age of the universe.
#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <cmath>
#include <cstdio>

using namespace SMVacuumDecay;

static int failures = 0;

// Locate the stability boundary in Mh at fixed Mt by bisection:
// for Mh below the boundary the vacuum is (meta)stable-negative... precisely:
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
    const double literature = 128.6; // Buttazzo et al. (2013) at Mt ~ 173.3 GeV
    std::printf("measured Mh_crit = %.3f GeV   literature ~ %.1f GeV\n", Mh_crit, literature);
    if (std::abs(Mh_crit - literature) > 2.5) { std::printf("FAIL boundary-band\n"); ++failures; }

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
