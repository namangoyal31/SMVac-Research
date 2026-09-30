// Emits convergence-study data as CSV on stdout, for the repository figures:
//   --mode dt      : optimized ansatz action vs RK4 step at the SM point
//   --mode quad    : action at the optimal radius vs Simpson points N
//
// The RG-trajectory construction intentionally replicates the inner loop of
// classify_stability (about 15 lines) so that this tool stays a leaf of the
// dependency graph; tests/physics/test_convergence.cpp enforces the same
// settings against the production pipeline.
//
// Usage: convergence_scan --mode dt|quad [--mh 125.1] [--mt 173.1]

#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <SMVacuumDecay/EffectivePotential.hpp>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>

using namespace SMVacuumDecay;

namespace {

void build_table(RGEHelper& rge, double Mh, double Mt, double dt) {
    rge.clear();
    StandardModelParameters y = get_nnlo_matching(Mh, Mt);
    double t = 2.0 * std::log(Mt);
    double tPlanck = 2.0 * std::log(planck_mass);
    rge.add_point(t, y);
    while (t < tPlanck) {
        if (t + dt > tPlanck) dt = tPlanck - t;
        y = rk4_single_step(y, t, dt);
        t += dt;
        rge.add_point(t, y);
    }
}

// Replicates the stability scan of classify_stability and returns the
// quantities needed for the action evaluation.
void locate_benchmark_scales(RGEHelper& rge, double Mh, double Mt,
                             double& mu_inst, double& t_min_lambda) {
    StandardModelParameters y = get_nnlo_matching(Mh, Mt);
    double t = 2.0 * std::log(Mt);
    double tPlanck = 2.0 * std::log(planck_mass);
    double dt = 0.1;
    double t_v = 2.0 * std::log(v);
    double min_lam_eff = 1.0;
    mu_inst = planck_mass;
    t_min_lambda = t;
    rge.add_point(t, y);
    while (t < tPlanck) {
        if (t + dt > tPlanck) dt = tPlanck - t;
        y = rk4_single_step(y, t, dt);
        t += dt;
        rge.add_point(t, y);
        if (t >= t_v) {
            double lam_eff = get_lambda_eff(y);
            if (lam_eff < min_lam_eff) { min_lam_eff = lam_eff; t_min_lambda = t; }
            if (lam_eff < 0 && mu_inst == planck_mass) mu_inst = std::exp(t / 2.0);
        }
    }
}

} // namespace

int main(int argc, char* argv[]) {
    std::string mode = "dt";
    double Mh = 125.1, Mt = 173.1;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--mode" && i + 1 < argc) mode = argv[++i];
        else if (arg == "--mh" && i + 1 < argc) Mh = std::stod(argv[++i]);
        else if (arg == "--mt" && i + 1 < argc) Mt = std::stod(argv[++i]);
        else { std::cerr << "Unknown argument: " << arg << std::endl; return 1; }
    }

    std::cout << std::setprecision(12);

    if (mode == "dt") {
        std::cout << "dt,S_exact\n";
        for (double dt : {0.4, 0.2, 0.1, 0.05, 0.025, 0.0125}) {
            StabilityResult r = classify_stability(Mh, Mt, dt);
            std::cout << dt << "," << r.S_exact << "\n";
        }
        return 0;
    }

    if (mode == "quad") {
        RGEHelper rge;
        double mu_inst, t_min_lambda;
        locate_benchmark_scales(rge, Mh, Mt, mu_inst, t_min_lambda);
        double R_opt = -1.0;
        find_minimum_action(rge, mu_inst, t_min_lambda, &R_opt);
        std::cout << "N,S_at_Ropt\n";
        for (int N : {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384}) {
            ActionEvaluation e = evaluate_action_components(rge, mu_inst, R_opt, N);
            std::cout << N << "," << e.total << "\n";
        }
        return 0;
    }

    std::cerr << "Unknown mode: " << mode << std::endl;
    return 1;
}
