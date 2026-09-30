#include <SMVacuumDecay/EffectivePotential.hpp>
#include <cmath>

namespace SMVacuumDecay {

double get_lambda_eff(const StandardModelParameters& p) {
    // One-loop Coleman-Weinberg correction to the quartic coupling, evaluated
    // with running couplings at mu = phi, so that the logs are the standard
    // field-dependent mass ratios m_i^2(phi)/phi^2:
    //   m_t^2 = yt^2 phi^2 / 2,  m_W^2 = g2^2 phi^2 / 4,
    //   m_Z^2 = (g1^2 + g2^2) phi^2 / 4.
    // In terms of V1 = sum_i n_i m_i^4/(64 pi^2) (ln(m_i^2/mu^2) - c_i) with
    // mu = phi, lambda_eff = lambda + 4 V1 / phi^4.
    double yt2 = p.yt * p.yt;
    double yt4 = yt2 * yt2;
    double term_t = -3.0 * yt4 * (std::log(0.5 * yt2) - 1.5);

    double g2_2 = p.g2 * p.g2;
    double g2_4 = g2_2 * g2_2;
    double term_W = 0.375 * g2_4 * (std::log(0.25 * g2_2) - 5.0/6.0);

    double g12 = p.g1 * p.g1 + g2_2;
    double g12_2 = g12 * g12;
    double term_Z = 0.1875 * g12_2 * (std::log(0.25 * g12) - 5.0/6.0);

    double delta_lambda = (term_t + term_W + term_Z) / (16.0 * PI2);
    return p.lambda + delta_lambda;
}

double V_eff(double phi, RGEHelper& rge) {
    if (phi <= 0) return 0.0;
    // mu = phi (RG improvement at the field scale); t = ln(mu^2).
    double t = 2.0 * std::log(phi);
    StandardModelParameters p = rge.get_params(t);
    double lam = get_lambda_eff(p);
    return 0.25 * lam * phi * phi * phi * phi;
}

}
