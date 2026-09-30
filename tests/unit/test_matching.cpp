// Literature check of the NNLO matching conditions: the implemented
// linearized matching is evaluated at the central masses and compared with
// the published values of Buttazzo et al. (2013), eqs. (55)-(60) and
// Table 3: lambda(Mt) = 0.12604, yt(Mt) = 0.93690, g2(Mt) = 0.64779,
// g3(Mt) = 1.1666, g1(Mt) = sqrt(5/3) g_Y = 0.46246.
#include <SMVacuumDecay/RGE.hpp>
#include <cmath>
#include <cstdio>

using namespace SMVacuumDecay;

static int failures = 0;

#define CHECK_CLOSE(val, ref, tol, name) do { \
    double d_ = std::abs((val) - (ref)); \
    std::printf("%-28s computed %.8f  published %.8f  diff %.2e\n", name, (double)(val), (double)(ref), d_); \
    if (!(d_ < (tol))) { std::printf("FAIL %s\n", name); ++failures; } \
} while (0)

int main() {
    // Central masses assumed by the published linearized matching.
    double Mh = 125.15, Mt = 173.34;
    StandardModelParameters p = get_nnlo_matching(Mh, Mt);

    CHECK_CLOSE(p.lambda, 0.12604, 1e-4, "lambda(Mt) vs Buttazzo");
    CHECK_CLOSE(p.yt,     0.93690, 1e-4, "yt(Mt) vs Buttazzo");
    CHECK_CLOSE(p.g2,     0.64779, 1e-9, "g2(Mt) vs Buttazzo");
    CHECK_CLOSE(p.g3,     1.1666,  1e-9, "g3(Mt) vs Buttazzo");
    CHECK_CLOSE(p.g1,     std::sqrt(5.0/3.0) * 0.35830, 1e-6, "g1(Mt) = sqrt(5/3) g_Y");

    // Tree-level consistency: the matching reproduces Mh^2 = 2 lambda v^2 up
    // to the NNLO offset (-0.00313), by construction of the parameterization.
    double lambda_tree = Mh * Mh / (2.0 * v * v);
    CHECK_CLOSE(p.lambda - lambda_tree, -0.00313, 1e-12, "NNLO offset");

    if (failures == 0) { std::printf("PASS\n"); return 0; }
    std::printf("FAIL (%d)\n", failures);
    return 1;
}
