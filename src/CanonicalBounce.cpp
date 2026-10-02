#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <SMVacuumDecay/EffectivePotential.hpp>
#include <SMVacuumDecay/numerics/GoldenSectionSearch.hpp>
#include <SMVacuumDecay/numerics/Simpson.hpp>
#include <cmath>
#include <algorithm>
#include <limits>

namespace SMVacuumDecay {

namespace {
// Steepness of the compactified radial variable r = R e^{alpha x}, u = tanh x.
// The profile is exactly the Fubini-Lipatov form for any alpha; alpha = 4
// concentrates the quadrature points near the bubble wall.
constexpr double kUmapAlpha = 4.0;
// Number of Simpson subintervals for the potential integral over u in (-1, 1).
constexpr int kQuadraturePoints = 2048;
// For |lambda_R| <~ 3.3e-4 the kinetic term alone exceeds this bound; the
// potential integral is skipped and the kinetic term returned as an
// approximation. Such points exceed the metastability threshold by orders of
// magnitude either way (docs/audit-notes.md, finding 8); the returned points
// are flagged MethodKineticShortcut.
constexpr double kKineticShortcutThreshold = 5e5;
// Search bracket and tolerance of the radius optimization: +-8 in ln R,
// i.e. ~7 decades in R.
constexpr double kLogRadiusBracket = 8.0;
constexpr double kLogRadiusTolerance = 1e-13;
// Sentinel returned by the internal action evaluation when no radius in the
// search bracket gives lambda_R < 0. Used only inside the golden-section
// search; classify_stability converts it to MethodFence with S_trial = NaN
// before anything is exported.
constexpr double kFenceSentinel = 1e100;
} // namespace

double get_pure_sm_lambda(RGEHelper& rge, double phi_dimless, double mu_inst) {
    double h = phi_dimless * mu_inst;
    // Couplings are frozen above the Planck mass: the RGEs are not trusted
    // beyond it and the trial action is regularized there by the phi^6 term.
    double h_rge = std::min(h, planck_mass);
    if (h_rge <= v) return 0.0;

    double t = 2.0 * std::log(h_rge);
    StandardModelParameters p = rge.get_params(t);
    return get_lambda_eff(p);
}

double integrand_u(RGEHelper& rge, double u, double R, double mu_inst,
                   double prefactor, double c6) {
    if (std::abs(u) >= 1.0) return 0.0;
    double ratio = (1.0 + u) / (1.0 - u);
    double e2x   = std::pow(ratio, kUmapAlpha);   // = (r/R)^2
    double e4x   = e2x * e2x;                     // = (r/R)^4
    double phi_x  = prefactor * 2.0 / (R * (e2x + 1.0));

    double lam_phi = get_pure_sm_lambda(rge, phi_x, mu_inst);
    double phi4 = std::pow(phi_x, 4);

    double V_SM = 0.25 * lam_phi * phi4;
    // Planck-suppressed dimension-6 operator with free coefficient c6
    // (default kDefaultC6; an assumed EFT regulator, not a prediction --
    // docs/limitations.md, item 3).
    double V_6 = (c6 / 6.0) * std::pow(phi_x * mu_inst / planck_mass, 2) * phi4;

    double V_x = V_SM + V_6;

    // dr/du divided by r^4: with r = R e^{alpha x} and u = tanh x,
    // dr = r * alpha dx = r * alpha du/(1 - u^2), so
    // 2 pi^2 r^3 V dr = 2 pi^2 R^4 (r/R)^4 V * alpha/(1-u^2) du.
    double jacobian = kUmapAlpha / (1.0 - u * u);
    return 2.0 * pi * pi * std::pow(R, 4) * e4x * V_x * jacobian;
}

ActionEvaluation evaluate_action_components(RGEHelper& rge, double mu_inst, double R,
                                            int quad_points, double c6) {
    if (R <= 0.0) return {kFenceSentinel, kFenceSentinel, 0.0};
    double mu_R = mu_inst / R;
    if (mu_R <= 0.0) return {kFenceSentinel, kFenceSentinel, 0.0};
    double t_R    = 2.0 * std::log(mu_R);
    StandardModelParameters p_R    = rge.get_params(t_R);
    double lambda_R = p_R.lambda;
    if (lambda_R >= 0.0) return {kFenceSentinel, kFenceSentinel, 0.0};
    double abs_lam  = std::abs(lambda_R);
    double prefactor = std::sqrt(2.0 / abs_lam);
    // For the conformal profile phi = A * 2R/(R^2 + r^2) with A^2 = 2/|lambda_R|
    // the kinetic integral is exactly (8/3) pi^2 A^2 = 16 pi^2/(3 |lambda_R|);
    // this holds for ANY potential because the profile is held fixed and only
    // its amplitude and scale are set by lambda_R.
    double kinetic_term = (16.0 * pi * pi) / (3.0 * abs_lam);
    if (kinetic_term > kKineticShortcutThreshold) return {kinetic_term, kinetic_term, 0.0};
    auto f = [&](double u) {
        return integrand_u(rge, u, R, mu_inst, prefactor, c6);
    };
    double potential_integral = numerics::simpson_integrate(-1.0, 1.0, quad_points, f);

    return {kinetic_term + potential_integral, kinetic_term, potential_integral};
}

double evaluate_action_at_R(RGEHelper& rge, double mu_inst, double R,
                            int quad_points, double c6) {
    return evaluate_action_components(rge, mu_inst, R, quad_points, c6).total;
}

double find_minimum_action(RGEHelper& rge, double mu_inst, double t_min_lambda,
                           double* R_opt_out, double c6) {
    // Initial bracket: the conformal estimate places the optimal radius at
    // R ~ mu_inst / mu(lambda_min); the bracket spans +-8 in ln R (~7 decades
    // in R) around it. The action is treated as unimodal on this bracket (a
    // heuristic; see docs/limitations.md).
    double mu_opt = std::exp(t_min_lambda / 2.0);
    double R_opt = mu_inst / mu_opt;
    double logR_opt = std::log(R_opt);
    double a = logR_opt - kLogRadiusBracket;
    double b = logR_opt + kLogRadiusBracket;
    double tol = kLogRadiusTolerance;

    auto f = [&](double logR) {
        return evaluate_action_at_R(rge, mu_inst, std::exp(logR),
                                    kDefaultQuadraturePoints, c6);
    };

    double lo, hi;
    double S_min = numerics::golden_section_search(a, b, tol, f, &lo, &hi);
    if (R_opt_out) *R_opt_out = std::exp(0.5 * (lo + hi));
    return S_min;
}

StabilityResult classify_stability(double Mh, double Mt, double rge_step, double c6) {
    const double nan_v = std::numeric_limits<double>::quiet_NaN();
    StabilityResult out{};
    out.Mh = Mh;
    out.status = 0; out.method_flag = MethodOK;
    out.S_trial = nan_v; out.S_conformal = nan_v;
    out.S_kinetic = nan_v; out.S_potential = nan_v; out.S_threshold = nan_v;
    out.mu_inst = -1.0; out.lambda_min = 1.0;

    StandardModelParameters y = get_nnlo_matching(Mh, Mt);

    double t0 = 2*std::log(Mt);
    double tPlanck = 2*std::log(planck_mass);

    RGEHelper rge;
    bool is_unstable = false;
    double mu_inst = planck_mass;

    double min_lambda_eff = 1.0;
    double t_min_lambda = t0;

    double t = t0;
    double dt = rge_step;
    rge.add_point(t, y);

    double t_v = 2 * std::log(v);

    if (t >= t_v) {
        min_lambda_eff = get_lambda_eff(y);
        if (min_lambda_eff <= 0.0) {
            mu_inst = std::exp(t/2.0);
            is_unstable = true;
        }
    }

    while (t < tPlanck) {
        if (t + dt > tPlanck) dt = tPlanck - t;

        y = rk4_single_step(y, t, dt);
        t += dt;
        rge.add_point(t, y);

        if (std::abs(y.lambda) > 4*pi || y.yt > 4*pi) {
            out.status = 4; out.method_flag = MethodPerturbativityLost;
            out.lambda_min = nan_v;
            return out;
        }
        if (!std::isfinite(y.lambda) || !std::isfinite(y.yt)) {
            out.status = 4; out.method_flag = MethodPerturbativityLost;
            out.lambda_min = nan_v;
            return out;
        }

        if (t >= t_v) {
            double lam_eff = get_lambda_eff(y);
            if (lam_eff < min_lambda_eff) {
                min_lambda_eff = lam_eff;
                t_min_lambda = t;
            }
            if (!is_unstable && lam_eff < 0) {
                mu_inst = std::exp(t/2.0);
                is_unstable = true;
            }
        }
    }

    out.lambda_min = min_lambda_eff;

    if (!is_unstable) {
        out.status = 1;  // Stable: lambda_eff never crosses zero
        out.method_flag = MethodOK;
        return out;
    }

    out.mu_inst = mu_inst;
    out.S_conformal = 8.0 * pi * pi / (3.0 * std::abs(min_lambda_eff));

    // Metastability threshold from the age of the universe:
    // Gamma/V ~ mu^4 e^{-S} <~ t_U^-4  <=>  S > 4 ln(mu * t_U),
    // with mu = the scale at which lambda_eff first crosses zero. This is a
    // conventional choice of the prefactor scale, not a determination of the
    // prefactor (docs/theory.md, docs/limitations.md).
    out.S_threshold = 4.0 * std::log(universe_age_GeV_minus1 * mu_inst);

    double R_opt = -1.0;
    double S_found = find_minimum_action(rge, mu_inst, t_min_lambda, &R_opt, c6);

    if (S_found >= 1e99) {
        // No radius in the search bracket gives lambda_R < 0: the trial
        // family produces no action here. The point is left UNCLASSIFIED.
        out.status = 0;
        out.method_flag = MethodFence;
        // S_trial / kinetic / potential stay NaN.
        return out;
    }

    ActionEvaluation best = evaluate_action_components(rge, mu_inst, R_opt,
                                                       kQuadraturePoints, c6);
    out.S_trial = best.total;
    out.S_kinetic = best.kinetic;
    out.S_potential = best.potential;

    if (out.S_trial <= 0.0) {
        // Non-positive trial action: the fixed-profile ansatz breaks down
        // (its amplitude sqrt(2/|lambda_R|) diverges as lambda_R -> 0^- while
        // the potential integrand uses lambda_eff). Not a physical verdict.
        out.status = 0;
        out.method_flag = MethodAnsatzFailed;
        return out;
    }

    out.method_flag = (best.kinetic >= kKineticShortcutThreshold)
                          ? MethodKineticShortcut : MethodOK;

    if (out.S_trial > out.S_threshold) { out.status = 2; return out; }  // Metastable
    out.status = 3;  // Unstable
    return out;
}

}
