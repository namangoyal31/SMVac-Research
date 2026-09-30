// Regression test: the RG-improved ansatz action at the low-Higgs point
// (Mh = 5 GeV, Mt = 105 GeV) must reproduce the frozen reference value.
// This point is deliberately pathological: far outside the conformal regime
// the ansatz action turns NEGATIVE (the potential integral overwhelms the
// kinetic term), which is a documented breakdown of the ansatz, not a
// numerical error (see docs/limitations.md). Regression, not independent
// physics validation.
#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <SMVacuumDecay/EffectivePotential.hpp>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace SMVacuumDecay;

int main() {
    double Low_Mh = 5.0;
    double Low_Mt = 105.0;
    double Ref_Numerical_Action = -1526742.427981734;

    StandardModelParameters y = get_nnlo_matching(Low_Mh, Low_Mt);

    double t0 = 2*std::log(Low_Mt);
    double tPlanck = 2*std::log(planck_mass);

    RGEHelper rge;
    bool is_unstable = false;
    double mu_inst = planck_mass;

    double min_lambda_eff = 1.0;
    double t_min_lambda = t0;

    double t = t0;
    double dt = 0.1;
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

        if (std::abs(y.lambda) > 4*pi || y.yt > 4*pi) break;
        if (!std::isfinite(y.lambda) || !std::isfinite(y.yt)) break;

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

    if (!is_unstable) {
        std::cout << "Test failed: Low-Mass point evaluated as stable!" << std::endl;
        return 1;
    }

    double S_exact = find_minimum_action(rge, mu_inst, t_min_lambda);

    double err = std::abs((S_exact - Ref_Numerical_Action) / Ref_Numerical_Action);

    std::cout << std::scientific << std::setprecision(12);
    std::cout << "Computed S_exact: " << S_exact << "\n";
    std::cout << "Reference S_exact: " << Ref_Numerical_Action << "\n";
    std::cout << "Bounce Regression Error: " << err << std::endl;

    if (err < 1e-12) {
        std::cout << "PASS" << std::endl;
        return 0;
    } else {
        std::cout << "FAIL" << std::endl;
        return 1;
    }
}
