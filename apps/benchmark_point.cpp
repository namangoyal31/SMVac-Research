// Single-point benchmark: evaluates both stability estimators at one
// (Mh, Mt) point and prints all diagnostic quantities.
//
// Usage: benchmark_point <Mh> <Mt>

#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <SMVacuumDecay/FubiniLipatov.hpp>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>

using namespace std;
using namespace SMVacuumDecay;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: benchmark_point <Mh> <Mt>" << endl;
        return 1;
    }

    double Mh = std::stod(argv[1]);
    double Mt = std::stod(argv[2]);

    cout << std::scientific << std::setprecision(12);
    cout << "--- Benchmark Point ---" << endl;
    cout << "Mh: " << Mh << " GeV" << endl;
    cout << "Mt: " << Mt << " GeV" << endl;

    auto fl_res = classify_buttazzo(Mh, Mt);
    cout << "\n--- Analytical (Fubini-Lipatov, strict conformal) ---" << endl;
    cout << "Status: " << std::get<0>(fl_res) << endl;
    cout << "S_approx: " << std::get<1>(fl_res) << endl;

    StabilityResult res = classify_stability(Mh, Mt);
    cout << "\n--- Numerical (RG-improved Fubini-Lipatov ansatz) ---" << endl;
    cout << "Status: " << res.status << endl;
    cout << "S_exact: " << res.S_exact << endl;
    cout << "S_approx: " << res.S_approx << endl;
    cout << "S_kinetic: " << res.S_kinetic << endl;
    cout << "S_potential: " << res.S_potential << endl;
    cout << "S_threshold: " << res.S_threshold << endl;
    cout << "mu_inst: " << res.mu_inst << endl;
    cout << "lambda_min: " << res.lambda_min << endl;

    return 0;
}
