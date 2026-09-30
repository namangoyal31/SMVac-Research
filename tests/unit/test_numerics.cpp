// Unit tests for the generic numerics: RK4 convergence order, Simpson
// quadrature exactness and convergence order, golden-section search.
// These validate the algorithms in isolation, independently of the physics.
#include <SMVacuumDecay/numerics/RK4.hpp>
#include <SMVacuumDecay/numerics/Simpson.hpp>
#include <SMVacuumDecay/numerics/GoldenSectionSearch.hpp>
#include <cmath>
#include <cstdio>
#include <functional>

using namespace SMVacuumDecay::numerics;

static const double kPi = 3.14159265358979323846;
static int failures = 0;

#define CHECK_LT(a, b, name) do { \
    if (!((a) < (b))) { \
        std::printf("FAIL %s: %g not < %g\n", name, (double)(a), (double)(b)); \
        ++failures; \
    } \
} while (0)

int main() {
    // --- RK4 order on the actual RGE test system dlambda/dt = beta ---
    // Integrate dy/dt = -y (exponential decay) from 0 to 10; error must
    // decrease by ~16x per halving (4th order).
    auto decay = [](double y, double) { return -y; };
    auto integrate = [&](double dt) {
        double y = 1.0;
        for (double t = 0; t < 10.0; t += dt) {
            double h = (t + dt > 10.0) ? 10.0 - t : dt;
            y = rk4_single_step(y, t, h, decay);
        }
        return y;
    };
    const double exact = std::exp(-10.0);
    double e1 = std::abs(integrate(0.1) - exact);
    double e2 = std::abs(integrate(0.05) - exact);
    double order = std::log2(e1 / e2);
    std::printf("RK4: err(0.1)=%.3e err(0.05)=%.3e observed order=%.2f\n", e1, e2, order);
    CHECK_LT(order, 5.0, "rk4-order-too-fast");
    CHECK_LT(3.0, order, "rk4-order-too-slow");

    // --- Simpson: cubic polynomials are integrated exactly ---
    auto f3 = [](double x) { return x*x*x - 2.0*x*x + 0.5*x + 3.0; };
    auto F3 = [](double x) { return x*x*x*x/4.0 - 2.0*x*x*x/3.0 + 0.25*x*x + 3.0*x; };
    double a = -1.5, b = 2.25;
    double num = simpson_integrate(a, b, 512, f3);
    double ana = F3(b) - F3(a);
    std::printf("Simpson cubic: num=%.15f ana=%.15f diff=%.3e\n", num, ana, std::abs(num - ana));
    CHECK_LT(std::abs(num - ana), 1e-12, "simpson-cubic-exact");

    // --- Simpson convergence order on a smooth non-polynomial ---
    auto fs = [](double x) { return std::sin(x); };
    double exact_s = 2.0; // integral of sin over [0, pi]
    double q1 = std::abs(simpson_integrate(0.0, kPi, 128, fs) - exact_s);
    double q2 = std::abs(simpson_integrate(0.0, kPi, 256, fs) - exact_s);
    double order_s = std::log2(q1 / q2);
    std::printf("Simpson sin: err(128)=%.3e err(256)=%.3e observed order=%.2f\n", q1, q2, order_s);
    CHECK_LT(3.0, order_s, "simpson-order");

    // --- Golden-section search on a parabola with a known minimum ---
    const double true_min_x = 1.2345;
    auto parab = [&](double logR) { return (logR - true_min_x) * (logR - true_min_x) + 7.0; };
    double lo, hi;
    double best = golden_section_search(-5.0, 8.0, 1e-13, parab, &lo, &hi);
    double min_x = 0.5 * (lo + hi);
    std::printf("GSS: x=%.10f (true %.10f) f=%.12f\n", min_x, true_min_x, best);
    CHECK_LT(std::abs(min_x - true_min_x), 1e-8, "gss-location");
    CHECK_LT(std::abs(best - 7.0), 1e-12, "gss-value");

    if (failures == 0) { std::printf("PASS\n"); return 0; }
    std::printf("FAIL (%d)\n", failures);
    return 1;
}
