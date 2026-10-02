// Unit test: the Coleman-Weinberg Z-boson term must use the PHYSICAL
// hypercharge combination gY^2 + g2^2 = (3/5) g1^2 + g2^2, because
// StandardModelParameters stores the GUT-normalized g1 = sqrt(5/3) gY
// (m_Z^2 = (gY^2 + g2^2) phi^2 / 4). An earlier version used
// g1^2 + g2^2, overestimating gY^2 by 5/3 (docs/audit-notes.md).
//
// The test evaluates get_lambda_eff at a fixed parameter point and compares
// against both the correct and the incorrect closed-form expression, so it
// discriminates between the two normalizations rather than just mirroring
// the implementation.
#include <SMVacuumDecay/EffectivePotential.hpp>
#include <cmath>
#include <cstdio>

using namespace SMVacuumDecay;

static int failures = 0;

// Closed-form Coleman-Weinberg correction, with the Z mass combination
// passed explicitly so both normalizations can be compared.
static double closed_form(const StandardModelParameters& p, double zComb) {
    double yt2 = p.yt * p.yt;
    double term_t = -3.0 * yt2 * yt2 * (std::log(0.5 * yt2) - 1.5);
    double g2_2 = p.g2 * p.g2;
    double term_W = 0.375 * g2_2 * g2_2 * (std::log(0.25 * g2_2) - 5.0 / 6.0);
    double term_Z = 0.1875 * zComb * zComb * (std::log(0.25 * zComb) - 5.0 / 6.0);
    return p.lambda + (term_t + term_W + term_Z) / (16.0 * PI2);
}

int main() {
    // Representative running point (high-scale-like couplings).
    StandardModelParameters p{};
    p.g1 = 0.52;    // GUT-normalized (gY^2 = 0.6 * 0.52^2 = 0.1622)
    p.g2 = 0.55;
    p.g3 = 1.05;    // unused by lambda_eff
    p.yt = 0.42;
    p.yb = 0.02;
    p.ytau = 0.01;
    p.lambda = -0.01;

    double computed = get_lambda_eff(p);

    double zCorrect = 0.6 * p.g1 * p.g1 + p.g2 * p.g2;   // (3/5) g1^2 + g2^2
    double zWrong   = p.g1 * p.g1 + p.g2 * p.g2;         // GUT-norm. bug

    double expected_correct = closed_form(p, zCorrect);
    double expected_wrong   = closed_form(p, zWrong);

    std::printf("lambda_eff(computed)            = %.12f\n", computed);
    std::printf("lambda_eff(closed form, gY^2)   = %.12f\n", expected_correct);
    std::printf("lambda_eff(closed form, g1^2)   = %.12f  (the bug)\n", expected_wrong);
    std::printf("discrimination |correct - wrong| = %.3e\n",
                std::abs(expected_correct - expected_wrong));

    // The implementation must match the (3/5) g1^2 + g2^2 closed form...
    if (std::abs(computed - expected_correct) > 1e-14 * std::abs(expected_correct)) {
        std::printf("FAIL lambda-eff-correct-normalization\n"); ++failures;
    }
    // ...and must be clearly distinguishable from the buggy one.
    if (std::abs(expected_correct - expected_wrong) < 1e-6) {
        std::printf("FAIL test-insensitive (normalizations degenerate at this point)\n");
        ++failures;
    }
    if (std::abs(computed - expected_wrong) < 0.5 * std::abs(expected_correct - expected_wrong)) {
        std::printf("FAIL lambda-eff-still-using-g1^2\n"); ++failures;
    }

    // Structural check: with g1 -> 0 the Z and W field-dependent masses are
    // equal (m_Z^2 = g2^2 phi^2/4 = m_W^2), but the loop degeneracy factors
    // differ (n_W = 6, n_Z = 3), so the Z term is exactly HALF the W term.
    // This pins the degree-of-freedom counting of the Z contribution.
    StandardModelParameters p0 = p;
    p0.g1 = 0.0;
    double computed_g0 = get_lambda_eff(p0);
    double g2_2 = p.g2 * p.g2;
    double term_t0 = -3.0 * p.yt * p.yt * p.yt * p.yt * (std::log(0.5 * p.yt * p.yt) - 1.5);
    double term_w = 0.375 * g2_2 * g2_2 * (std::log(0.25 * g2_2) - 5.0 / 6.0);
    double term_z_half = 0.5 * term_w;  // n_Z/n_W = 1/2 at equal masses
    double expected_g0 = p.lambda + (term_w + term_z_half + term_t0) / (16.0 * PI2);
    std::printf("g1 = 0: computed %.12f vs closed form %.12f\n", computed_g0, expected_g0);
    if (std::abs(computed_g0 - expected_g0) > 1e-14 * std::abs(expected_g0)) {
        std::printf("FAIL g1-equals-zero-degeneracy\n"); ++failures;
    }

    if (failures == 0) { std::printf("PASS\n"); return 0; }
    std::printf("FAIL (%d)\n", failures);
    return 1;
}
