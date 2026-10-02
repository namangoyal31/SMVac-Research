// Physics test: the classification/method-flag separation and the c6
// sensitivity of the trial action.
//
// 1. Flag taxonomy: the SM benchmark point must be status 2 (metastable)
//    with MethodOK; the pathological low-Higgs point (5, 105) must be
//    status 0 (undetermined) with MethodAnsatzFailed and a negative trial
//    action -- never a physical "unstable" verdict; a fully stable coupling
//    table must drive the action evaluation into the fence branch.
// 2. c6 separation: at the near-boundary point (134.75, 176.5) the trial
//    action with the default c6 = 1 exceeds the conformal estimate by >20%,
//    while with c6 = 0 the two agree to <2%: the near-boundary action
//    difference is dominated by the assumed Planck-suppressed operator,
//    not by RG improvement (docs/audit-notes.md, referee report M-2).
#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <SMVacuumDecay/EffectivePotential.hpp>
#include <cmath>
#include <cstdio>

using namespace SMVacuumDecay;

static int failures = 0;

static void fill_constant_table(RGEHelper& rge, const StandardModelParameters& p) {
    rge.clear();
    for (double t = 0.0; t <= 2.0 * std::log(planck_mass); t += 0.05) {
        rge.add_point(t, p);
    }
}

int main() {
    std::printf("=== 1. Method flags ===\\n");
    StabilityResult sm = classify_stability(125.1, 173.1);
    std::printf("(125.1, 173.1): status %d flag %d S_trial = %.6f\\n",
                sm.status, sm.method_flag, sm.S_trial);
    if (sm.status != 2 || sm.method_flag != MethodOK) {
        std::printf("FAIL sm-ok\\n"); ++failures;
    }

    StabilityResult low = classify_stability(5.0, 105.0);
    std::printf("(5, 105): status %d flag %d S_trial = %.6e\\n",
                low.status, low.method_flag, low.S_trial);
    if (low.status != 0 || low.method_flag != MethodAnsatzFailed) {
        std::printf("FAIL low-ansatz-failed\\n"); ++failures;
    }
    if (!(low.S_trial < 0.0)) {
        std::printf("FAIL low-negative-action\\n"); ++failures;
    }

    // Fence branch: a constant POSITIVE lambda table leaves no radius with
    // lambda_R < 0; evaluate_action_components must return the internal
    // fence sentinel and find_minimum_action must propagate it.
    StandardModelParameters pstable{};
    pstable.g1 = 0.46; pstable.g2 = 0.65; pstable.g3 = 1.0;
    pstable.yt = 0.1; pstable.yb = 0.01; pstable.ytau = 0.01;
    pstable.lambda = +0.05;
    RGEHelper rge;
    fill_constant_table(rge, pstable);
    double fenceS = find_minimum_action(rge, 1e6, 2.0 * std::log(1e10));
    std::printf("stable-table find_minimum_action = %.3e (>= 1e100 expected)\\n", fenceS);
    if (!(fenceS >= 1e99)) { std::printf("FAIL fence-sentinel\\n"); ++failures; }

    std::printf("=== 2. c6 sensitivity at (134.75, 176.5) ===\\n");
    StabilityResult c0 = classify_stability(134.75, 176.5, kDefaultRgeStep, 0.0);
    StabilityResult c1 = classify_stability(134.75, 176.5);  // c6 = 1
    double conf = c1.S_conformal;
    double frac0 = c0.S_trial / conf - 1.0;
    double frac1 = c1.S_trial / conf - 1.0;
    std::printf("S_conformal = %.4f\\n", conf);
    std::printf("c6 = 0: S_trial = %.4f  (%+.2f%% vs conformal), flag %d\\n",
                c0.S_trial, 100 * frac0, c0.method_flag);
    std::printf("c6 = 1: S_trial = %.4f  (%+.2f%% vs conformal), flag %d\\n",
                c1.S_trial, 100 * frac1, c1.method_flag);
    // With c6 = 0 the RG improvement alone changes the action by < 2%...
    if (std::abs(frac0) > 0.02) { std::printf("FAIL c6-zero-rg-small\\n"); ++failures; }
    // ...while the assumed c6 = 1 operator dominates (> 20%).
    if (frac1 < 0.20) { std::printf("FAIL c6-one-dominates\\n"); ++failures; }
    // Both evaluations must be method-clean at this point.
    if (c0.method_flag != MethodOK || c1.method_flag != MethodOK) {
        std::printf("FAIL c6-points-method-ok\\n"); ++failures;
    }

    if (failures == 0) { std::printf("PASS\\n"); return 0; }
    std::printf("FAIL (%d)\\n", failures);
    return 1;
}
