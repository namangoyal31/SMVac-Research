// Regression test: both classifiers at the benchmark SM point must reproduce
// the frozen reference values of reference/v1.1/benchmark_results*.json.
// Regression, not independent physics validation (see tests/physics and
// docs/validation.md).
#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <SMVacuumDecay/FubiniLipatov.hpp>
#include <iostream>
#include <cmath>

using namespace SMVacuumDecay;

void assert_close(double actual, double expected, double tol, const char* name) {
    if (std::abs(actual - expected) > tol) {
        std::cerr << "FAIL: " << name << " mismatch. Actual: " << actual << " Expected: " << expected << std::endl;
        std::exit(1);
    }
}

void assert_equal(int actual, int expected, const char* name) {
    if (actual != expected) {
        std::cerr << "FAIL: " << name << " mismatch. Actual: " << actual << " Expected: " << expected << std::endl;
        std::exit(1);
    }
}

int main() {
    double Mh = 125.1;
    double Mt = 173.1;

    std::cout << "Running Regression Validation for Benchmark Point: Mh=" << Mh << " GeV, Mt=" << Mt << " GeV" << std::endl;

    // Fubini-Lipatov (Analytical) Benchmark
    auto fl_res = classify_buttazzo(Mh, Mt);
    int fl_status = std::get<0>(fl_res);
    double fl_action = std::get<1>(fl_res);

    // Canonical Bounce (Numerical) Benchmark
    StabilityResult num_res = classify_stability(Mh, Mt);

    // Expected values from reference/v1.1 (regenerated after the matching and
    // interpolation corrections; see docs/audit-notes.md and CHANGELOG.md)
    int expected_status = 2; // Metastable
    double expected_fl_action = 2051.1373669116429;
    double expected_num_action = 2103.1046353416655;

    std::cout << "\nValidating Analytical solver..." << std::endl;
    assert_equal(fl_status, expected_status, "FL Status");
    assert_close(fl_action, expected_fl_action, 1e-9, "FL Action");

    std::cout << "Validating Numerical solver..." << std::endl;
    assert_equal(num_res.status, expected_status, "Numerical Status");
    assert_close(num_res.S_exact, expected_num_action, 1e-9, "Numerical Exact Action");
    assert_close(num_res.S_approx, expected_fl_action, 1e-9, "Numerical Approx Action");

    std::cout << "\nAll regression validations PASS!" << std::endl;
    return 0;
}
