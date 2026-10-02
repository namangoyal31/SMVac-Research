// Single-point benchmark: evaluates both estimators at one (Mh, Mt) point
// and prints all diagnostic quantities, including the method flag.
//
// Usage: benchmark_point <Mh> <Mt> [--c6 value]

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
    double c6 = kDefaultC6;
    for (int i = 3; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--c6" && i + 1 < argc) c6 = std::stod(argv[++i]);
        else { std::cerr << "Unknown argument: " << arg << std::endl; return 1; }
    }

    cout << std::scientific << std::setprecision(12);
    cout << "--- Benchmark Point ---" << endl;
    cout << "Mh: " << Mh << " GeV" << endl;
    cout << "Mt: " << Mt << " GeV" << endl;

    auto fl_res = classify_conformal(Mh, Mt);
    cout << "\n--- Conformal (strict Fubini-Lipatov estimate) ---" << endl;
    cout << "Status: " << std::get<0>(fl_res) << endl;
    cout << "S_conformal: " << std::get<1>(fl_res) << endl;

    StabilityResult res = classify_stability(Mh, Mt, kDefaultRgeStep, c6);
    cout << "\n--- Trial profile (RG-improved Fubini-Lipatov ansatz, c6 = "
         << c6 << ") ---" << endl;
    cout << "Status: " << res.status << "  (method_flag: " << res.method_flag
         << "; 0=OK 1=shortcut 2=fence 3=ansatz-failed 4=pert-lost)" << endl;
    cout << "S_trial: " << res.S_trial << endl;
    cout << "S_conformal: " << res.S_conformal << endl;
    cout << "S_kinetic: " << res.S_kinetic << endl;
    cout << "S_potential: " << res.S_potential << endl;
    cout << "S_threshold: " << res.S_threshold << endl;
    cout << "mu_inst: " << res.mu_inst << endl;
    cout << "lambda_min: " << res.lambda_min << endl;

    return 0;
}
