// Physics validation against analytic limits.
//
// 1. Constant-coupling conformal limit: with a frozen coupling table
//    (lambda = const < 0), the Fubini-Lipatov profile IS the exact bounce
//    and the action has the closed form
//        S(R) = 16 pi^2/(3|lambda_R|) - 8 pi^2 |lambda_eff| / (3 lambda_R^2),
//    where lambda_R is the raw running coupling used for the amplitude and
//    lambda_eff = lambda + Delta_lambda_CW enters the potential integrand.
//    The code must reproduce this identity for every profile scale R, which
//    validates the profile, the compactified Jacobian, the Simpson
//    quadrature, and the closed-form kinetic term end-to-end.
// 2. Deep negative-coupling regime: at (Mh, Mt) = (115, 180) GeV the
//    running is slow over the bubble and the trial action must approach the
//    conformal estimate S_conformal = 8 pi^2/(3|lambda_min|) to better than
//    0.1%. (This point is classified UNSTABLE by the threshold criterion;
//    the test checks internal consistency of the two estimators, not a
//    physical label.)
// 3. Virial identity: in the conformal regime S_kinetic/|S_potential| = 2;
//    the RG improvement only slightly perturbs this at the SM point.
#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <SMVacuumDecay/EffectivePotential.hpp>
#include <cmath>
#include <cstdio>

using namespace SMVacuumDecay;

static int failures = 0;

#define CHECK_CLOSE(val, ref, rtol, name) do { \
    double d_ = std::abs((val) - (ref)) / std::abs(ref); \
    std::printf("%-42s computed %14.10f  expected %14.10f  rel.diff %.2e\n", name, (double)(val), (double)(ref), d_); \
    if (!(d_ < (rtol))) { std::printf("FAIL %s\n", name); ++failures; } \
} while (0)

// Build an RG table that is constant in t (frozen couplings), matching the
// production table layout.
static void fill_constant_table(RGEHelper& rge, const StandardModelParameters& p) {
    rge.clear();
    for (double t = 0.0; t <= 2.0 * std::log(planck_mass); t += 0.05) {
        rge.add_point(t, p);
    }
}

int main() {
    std::printf("=== 1. Constant-coupling conformal limit ===\n");
    // Small but non-zero Yukawas/gaugings keep the CW correction finite.
    StandardModelParameters p{};
    p.g1 = 0.1; p.g2 = 0.1; p.g3 = 1.0;
    p.yt = 0.0742; p.yb = 0.01; p.ytau = 0.01;
    p.lambda = -0.01;
    RGEHelper rge;
    fill_constant_table(rge, p);

    const double lambda_R = p.lambda;                      // raw coupling (amplitude)
    const double lambda_eff = get_lambda_eff(p);           // coupling in the potential
    const double S_kin_expected = 16.0 * PI2 / (3.0 * std::abs(lambda_R));
    const double S_tot_expected = S_kin_expected
        - 8.0 * PI2 * std::abs(lambda_eff) / (3.0 * lambda_R * lambda_R);

    std::printf("lambda_eff (frozen table) = %.12f (raw lambda = %.12f)\n", lambda_eff, lambda_R);

    const double mu_inst = 1e6;   // unit scale: phi(0) = 2 sqrt(2/|lambda_R|) mu_inst / R
    for (double R : {1e-3, 1e-2, 1e-1}) {
        ActionEvaluation e = evaluate_action_components(rge, mu_inst, R);
        CHECK_CLOSE(e.kinetic, S_kin_expected, 1e-12, "S_kinetic (closed form)");
        CHECK_CLOSE(e.total, S_tot_expected, 1e-8, "S(R) analytic identity");
    }
    // Scale invariance: the action must be R-independent in this limit.
    ActionEvaluation e1 = evaluate_action_components(rge, mu_inst, 1e-3);
    ActionEvaluation e2 = evaluate_action_components(rge, mu_inst, 1e-1);
    CHECK_CLOSE(e2.total, e1.total, 1e-8, "R-independence (scale invariance)");

    std::printf("=== 2. Deep negative-coupling regime: S_trial -> S_conformal ===\n");
    StabilityResult deep = classify_stability(115.0, 180.0);
    std::printf("(115,180): S_trial=%.6f S_conformal=%.6f status=%d\n",
                deep.S_trial, deep.S_conformal, deep.status);
    CHECK_CLOSE(deep.S_trial, deep.S_conformal, 1e-3, "S_trial vs S_conformal (deep)");

    std::printf("=== 3. Virial ratio at the SM point ===\n");
    StabilityResult sm = classify_stability(125.1, 173.1);
    double virial = sm.S_kinetic / std::abs(sm.S_potential);
    std::printf("(125.1,173.1): S_kinetic=%.6f S_potential=%.6f virial=%.6f\n",
                sm.S_kinetic, sm.S_potential, virial);
    // The virial ratio is exactly 2 in the pure-quartic limit. The running
    // couplings perturb it by a few percent -- a measure of the trial
    // profile's non-stationarity under the full potential (ansatz
    // systematics, docs/limitations.md). The tolerance below checks the
    // near-cancellation, not stationarity.
    if (std::abs(virial - 2.0) > 0.15) { std::printf("FAIL virial-ratio\n"); ++failures; }

    if (failures == 0) { std::printf("PASS\n"); return 0; }
    std::printf("FAIL (%d)\n", failures);
    return 1;
}
