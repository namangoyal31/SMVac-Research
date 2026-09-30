#include <SMVacuumDecay/FubiniLipatov.hpp>
#include <SMVacuumDecay/RGE.hpp>
#include <SMVacuumDecay/EffectivePotential.hpp>
#include <cmath>

namespace SMVacuumDecay {

// Analytical (strict conformal) estimate: the RG equations are run up to the
// Planck scale, the minimum lambda_eff is located, and the Fubini-Lipatov
// action S = 8 pi^2 / (3 |lambda_min|) is formed from it. This coincides with
// the exact bounce action in the pure-quartic limit, where the couplings do
// not run; away from it the RG-improved ansatz in CanonicalBounce.cpp is the
// more accurate estimator.
std::tuple<int, double, double> classify_buttazzo(double Mh_input, double Mt) {
    StandardModelParameters y = get_nnlo_matching(Mh_input, Mt);
    double Mh_calc = Mh_input;
    double t0 = 2*std::log(Mt);
    // Standard analytic plots only evaluate up to Planck scale
    double tMax = 2*std::log(planck_mass);

    double min_lambda_eff = 1.0;
    double mu1 = -1.0;
    bool went_negative = false;
    double t = t0;
    // Fixed-step RK4 in t = ln(mu^2); the step error is quantified in
    // docs/validation.md.
    double dt = 0.1;

    double t_v = 2 * std::log(v);

    if (t >= t_v) {
        min_lambda_eff = get_lambda_eff(y);
        if (min_lambda_eff <= 0.0) {
            went_negative = true;
            mu1 = std::exp(t / 2.0);
        }
    }

    while (t < tMax) {
        if (t + dt > tMax) dt = tMax - t;
        y = rk4_single_step(y, t, dt);
        t += dt;

        if (t >= t_v) {
            double lam_eff = get_lambda_eff(y);
            if (!went_negative && lam_eff <= 0.0) {
                mu1 = std::exp(t / 2.0);
                went_negative = true;
            }
            if (lam_eff < min_lambda_eff) min_lambda_eff = lam_eff;
        }
        if (std::abs(y.lambda) > 4*pi || y.yt > 4*pi || y.g1 > 2.0) break;
    }

    if (t < 2*std::log(planck_mass)) return std::make_tuple(4, -1.0, Mh_calc);
    if (min_lambda_eff >= 0.0) return std::make_tuple(1, -1.0, Mh_calc);

    double S_approx = 8.0 * pi * pi / (3.0 * std::abs(min_lambda_eff));
    double S_threshold = 4.0 * std::log(universe_age_GeV_minus1 * mu1);

    bool is_metastable_S = (S_approx > S_threshold);
    if (is_metastable_S) return std::make_tuple(2, S_approx, Mh_calc);
    else return std::make_tuple(3, S_approx, Mh_calc);
}

}
