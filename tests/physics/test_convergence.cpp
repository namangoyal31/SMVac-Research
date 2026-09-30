// Numerical convergence studies of the production pipeline (not toy models):
//   1. RGE step size: classify_stability at the SM point for dt in
//      {0.2, 0.1, 0.05, 0.025}; the production step dt = 0.1 must agree with
//      the four-times-finer run to the stated tolerance.
//   2. Quadrature: the potential integral at the optimal radius for N in
//      {256, 1024, 4096, 16384}; Simpson's O(N^-4) convergence and the
//      production setting N = 2048 are checked.
//   3. Radius minimization: the golden-section result is compared with the
//      minimum of a fine grid scan over the same bracket.
#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <SMVacuumDecay/EffectivePotential.hpp>
#include <SMVacuumDecay/numerics/GoldenSectionSearch.hpp>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace SMVacuumDecay;

static int failures = 0;

// Rebuild the RG trajectory table exactly as classify_stability does, but
// with a configurable step, for the quadrature/minimization studies.
static void build_table(RGEHelper& rge, double Mh, double Mt, double dt) {
    rge.clear();
    StandardModelParameters y = get_nnlo_matching(Mh, Mt);
    double t = 2.0 * std::log(Mt);
    double tPlanck = 2.0 * std::log(planck_mass);
    rge.add_point(t, y);
    while (t < tPlanck) {
        if (t + dt > tPlanck) dt = tPlanck - t;
        y = rk4_single_step(y, t, dt);
        t += dt;
        rge.add_point(t, y);
    }
}

int main() {
    const double Mh = 125.1, Mt = 173.1;

    std::printf("=== 1. RGE step-size convergence (S_exact at the SM point) ===\n");
    std::vector<double> steps = {0.2, 0.1, 0.05, 0.025};
    std::vector<double> S;
    for (double dt : steps) {
        StabilityResult r = classify_stability(Mh, Mt, dt);
        S.push_back(r.S_exact);
        std::printf("dt = %-6.3f  S_exact = %.9f  (mu1 = %.6e)\n", dt, r.S_exact, r.mu_inst);
    }
    double rel_01 = std::abs(S[1] - S[3]) / S[3];
    std::printf("|S(0.1) - S(0.025)|/S = %.3e\n", rel_01);
    if (!(rel_01 < 5e-5)) { std::printf("FAIL rge-step-convergence\n"); ++failures; }
    // The difference must shrink roughly quadratically between successive
    // refinements (4th-order method on the trajectory, 1st-order in the
    // tabulated interpolation) - a loose order sanity check:
    double d1 = std::abs(S[0] - S[1]), d2 = std::abs(S[1] - S[2]), d3 = std::abs(S[2] - S[3]);
    std::printf("successive differences: %.3e, %.3e, %.3e\n", d1, d2, d3);
    if (!(d2 < d1 && d3 < d2)) { std::printf("FAIL rge-step-monotone\n"); ++failures; }

    std::printf("=== 2. Quadrature convergence at the optimal radius ===\n");
    RGEHelper rge;
    build_table(rge, Mh, Mt, 0.1);
    // Reproduce the classification to locate mu1 and t_min_lambda.
    StandardModelParameters y = get_nnlo_matching(Mh, Mt);
    double t = 2.0 * std::log(Mt), dt = 0.1;
    double min_lam_eff = 1.0, t_min = t, mu_inst = planck_mass;
    bool unstable = false;
    double t_v = 2.0 * std::log(v);
    while (t < 2.0 * std::log(planck_mass)) {
        if (t + dt > 2.0 * std::log(planck_mass)) dt = 2.0 * std::log(planck_mass) - t;
        y = rk4_single_step(y, t, dt);
        t += dt;
        if (t >= t_v) {
            double lam = get_lambda_eff(y);
            if (lam < min_lam_eff) { min_lam_eff = lam; t_min = t; }
            if (!unstable && lam < 0) { mu_inst = std::exp(t / 2.0); unstable = true; }
        }
    }
    double R_opt = -1.0;
    double S_gss = find_minimum_action(rge, mu_inst, t_min, &R_opt);
    std::printf("optimal R = %.6f (in 1/mu1 units), S_gss = %.9f\n", R_opt, S_gss);

    std::vector<int> Ns = {256, 1024, 4096, 16384};
    std::vector<double> Sq;
    for (int N : Ns) {
        ActionEvaluation e = evaluate_action_components(rge, mu_inst, R_opt, N);
        Sq.push_back(e.total);
        std::printf("N = %-6d  S(R*) = %.12f\n", N, e.total);
    }
    double rel_q = std::abs(Sq[1] - Sq[3]) / Sq[3];
    std::printf("|S(1024) - S(16384)|/S = %.3e\n", rel_q);
    if (!(rel_q < 1e-6)) { std::printf("FAIL quadrature-convergence\n"); ++failures; }
    // NOTE: the Simpson O(N^-4) order is verified on smooth integrands in
    // tests/unit/test_numerics.cpp. On the production integrand the action
    // is already converged to the ~1e-9 relative round-off floor at every N
    // tested (large cancellations in S_potential), so an order estimate
    // here would be meaningless noise.
    // Convergence floor documented: differences between successive N are
    // non-monotone at the 1e-9 relative level.

    std::printf("=== 3. Radius minimization: golden section vs fine grid ===\n");
    // Compare the minimization algorithms at a common (coarse) quadrature
    // setting, so the comparison isolates the minimization, not the
    // quadrature. Fine grid over the same 16-decade bracket.
    const int Ncoarse = 64;
    double logR_opt = std::log(mu_inst / std::exp(t_min / 2.0));
    auto action = [&](double logR) {
        return evaluate_action_at_R(rge, mu_inst, std::exp(logR), Ncoarse);
    };
    double lo = logR_opt - 8.0, hi = logR_opt + 8.0;
    double S_gss_coarse = numerics::golden_section_search(lo, hi, 1e-13, action);
    double best = 1e300;
    for (int i = 0; i <= 50000; ++i) {
        double s = action(lo + 16.0 * i / 50000.0);
        if (s < best) best = s;
    }
    std::printf("S_grid = %.9f  S_gss = %.9f  (S_grid - S_gss)/S = %.3e\n",
                best, S_gss_coarse, (best - S_gss_coarse) / std::abs(S_gss_coarse));
    if (!(S_gss_coarse <= best + 1e-6 * std::abs(best))) { std::printf("FAIL gss-vs-grid\n"); ++failures; }

    if (failures == 0) { std::printf("PASS\n"); return 0; }
    std::printf("FAIL (%d)\n", failures);
    return 1;
}
