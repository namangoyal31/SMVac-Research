#pragma once
#include <cmath>

namespace SMVacuumDecay {
namespace numerics {

// Golden-section minimization of f over [a, b], assuming f is unimodal there.
// If a_out/b_out are provided they receive the final bracket, whose midpoint
// approximates the minimizer to the requested tolerance.
template <typename Func>
double golden_section_search(double a, double b, double tol, Func f, double* a_out = nullptr, double* b_out = nullptr) {
    const double invphi = (std::sqrt(5.0) - 1.0) / 2.0;
    const double invphi2 = (3.0 - std::sqrt(5.0)) / 2.0;

    double h = b - a;
    if (h <= tol) {
        if (a_out) *a_out = a;
        if (b_out) *b_out = b;
        return f(0.5 * (a + b));
    }

    int n = static_cast<int>(std::ceil(std::log(tol / h) / std::log(invphi)));

    double c = a + invphi2 * h;
    double d = a + invphi * h;

    double yc = f(c);
    double yd = f(d);

    for (int k = 0; k < n; ++k) {
        if (yc < yd) {
            b = d;
            d = c;
            yd = yc;
            h = invphi * h;
            c = a + invphi2 * h;
            yc = f(c);
        } else {
            a = c;
            c = d;
            yc = yd;
            h = invphi * h;
            d = a + invphi * h;
            yd = f(d);
        }
    }

    if (a_out) *a_out = a;
    if (b_out) *b_out = b;
    return f(0.5 * (a + b));
}

} // namespace numerics
} // namespace SMVacuumDecay
