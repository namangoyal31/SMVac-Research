// Unit test of the metastability-threshold constants and formula:
// Gamma/V ~ mu^4 e^{-S} <~ t_U^{-4}  <=>  S > 4 ln(mu t_U),
// with t_U = 10 Gyr expressed in GeV^-1. The expected value is recomputed
// here from first principles (independent of the header constants).
#include <SMVacuumDecay/RGE.hpp>
#include <cmath>
#include <cstdio>

using namespace SMVacuumDecay;

static int failures = 0;

#define CHECK_CLOSE(val, ref, rtol, name) do { \
    double d_ = std::abs((val) - (ref)) / std::abs(ref); \
    std::printf("%-34s computed %.10e  expected %.10e  rel.diff %.2e\n", name, (double)(val), (double)(ref), d_); \
    if (!(d_ < (rtol))) { std::printf("FAIL %s\n", name); ++failures; } \
} while (0)

int main() {
    // Independent derivation: 10 Gyr = 10^10 * 365.25 d (using 365.25-day
    // years here on purpose: the header uses 365 d, so allow 1e-3 relative).
    double t_U_independent = 10.0e9 * 365.25 * 86400.0 / 6.582e-25; // GeV^-1
    CHECK_CLOSE(universe_age_GeV_minus1, t_U_independent, 1e-3, "universe age [GeV^-1]");

    // Threshold formula at a representative zero-crossing scale mu1 = 10^10 GeV:
    double mu1 = 1e10;
    double S_th = 4.0 * std::log(universe_age_GeV_minus1 * mu1);
    double S_th_expected = 4.0 * std::log(4.7913e41 * 1e10);
    CHECK_CLOSE(S_th, S_th_expected, 1e-3, "S_threshold at mu1 = 1e10 GeV");

    // Sanity: the threshold must grow logarithmically with mu and exceed the
    // SM-relevant range (~300-600) for mu1 in [1e8, 1e13] GeV.
    double lo = 4.0 * std::log(universe_age_GeV_minus1 * 1e8);
    double hi = 4.0 * std::log(universe_age_GeV_minus1 * 1e13);
    std::printf("threshold range for mu1 in [1e8, 1e13] GeV: [%.1f, %.1f]\n", lo, hi);
    if (!(lo < hi && lo > 300.0 && hi < 650.0)) { std::printf("FAIL threshold-range\n"); ++failures; }

    // Monotonicity in mu.
    if (!(4.0 * std::log(universe_age_GeV_minus1 * 2e10) > S_th)) {
        std::printf("FAIL threshold-monotonic\n"); ++failures;
    }

    if (failures == 0) { std::printf("PASS\n"); return 0; }
    std::printf("FAIL (%d)\n", failures);
    return 1;
}
