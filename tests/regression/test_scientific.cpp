// Regression test: both estimators at the benchmark SM point must reproduce
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

    // Strict conformal estimate
    auto fl_res = classify_conformal(Mh, Mt);
    int fl_status = std::get<0>(fl_res);
    double fl_action = std::get<1>(fl_res);

    // RG-improved trial-profile action
    StabilityResult num_res = classify_stability(Mh, Mt);

    // Expected values from reference/v1.1 (regenerated after the U(1)
    // normalization fix, the 4-loop sign fix and the method-flag refactor;
    // see docs/audit-notes.md and CHANGELOG.md)
    int expected_status = 2; // Metastable
    double expected_conformal = 2120.340020693041;
    double expected_trial = 2168.6895796842282;

    std::cout << "Validating conformal estimator..." << std::endl;
    assert_equal(fl_status, expected_status, "Conformal Status");
    assert_close(fl_action, expected_conformal, 1e-9, "Conformal Action");

    std::cout << "Validating trial-profile estimator..." << std::endl;
    assert_equal(num_res.status, expected_status, "Trial Status");
    assert_equal(num_res.method_flag, MethodOK, "Trial MethodFlag");
    assert_close(num_res.S_trial, expected_trial, 1e-9, "Trial Action");
    assert_close(num_res.S_conformal, expected_conformal, 1e-9, "Trial-path Conformal Action");

    std::cout << "All regression validations PASS!" << std::endl;
    return 0;
}
