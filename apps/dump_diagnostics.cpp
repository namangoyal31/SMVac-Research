// Dumps diagnostic CSV files for the repository figures at a single
// (Mh, Mt) point:
//   trajectory.csv   : scale_mu, lambda_raw, lambda_eff along the RG trajectory
//   action_curve.csv : R_over_mu1, S_total, S_kinetic, S_potential
//
// Usage: dump_diagnostics --mh 125.1 --mt 173.1 --out-dir results/diagnostics

#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <SMVacuumDecay/EffectivePotential.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#endif

using namespace SMVacuumDecay;

int main(int argc, char* argv[]) {
    double Mh = 125.1, Mt = 173.1;
    std::string out_dir = "results/diagnostics";
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--mh" && i + 1 < argc) Mh = std::stod(argv[++i]);
        else if (arg == "--mt" && i + 1 < argc) Mt = std::stod(argv[++i]);
        else if (arg == "--out-dir" && i + 1 < argc) out_dir = argv[++i];
        else { std::cerr << "Unknown argument: " << arg << std::endl; return 1; }
    }
#ifdef _WIN32
    _mkdir(out_dir.c_str());
#else
    mkdir(out_dir.c_str(), 0755);
#endif

    std::cout << std::setprecision(12);

    // --- RG trajectory with lambda_eff ---
    {
        std::ofstream out(out_dir + "/trajectory.csv");
        if (!out.is_open()) { std::cerr << "cannot open trajectory.csv" << std::endl; return 1; }
        out << std::setprecision(12);
        out << "scale_mu,lambda_raw,lambda_eff\n";
        StandardModelParameters y = get_nnlo_matching(Mh, Mt);
        double t = 2.0 * std::log(Mt);
        double tPlanck = 2.0 * std::log(planck_mass);
        double dt = 0.1;
        auto emit = [&](double t, const StandardModelParameters& p) {
            double mu = std::exp(t / 2.0);
            out << mu << "," << p.lambda << "," << get_lambda_eff(p) << "\n";
        };
        emit(t, y);
        while (t < tPlanck) {
            if (t + dt > tPlanck) dt = tPlanck - t;
            y = rk4_single_step(y, t, dt);
            t += dt;
            emit(t, y);
        }
    }

    // --- Action curve S(R) with its kinetic/potential decomposition ---
    {
        // Build the trajectory table (same settings as classify_stability).
        RGEHelper rge;
        {
            StandardModelParameters y = get_nnlo_matching(Mh, Mt);
            double t = 2.0 * std::log(Mt);
            double tPlanck = 2.0 * std::log(planck_mass);
            double dt = 0.1;
            rge.add_point(t, y);
            while (t < tPlanck) {
                if (t + dt > tPlanck) dt = tPlanck - t;
                y = rk4_single_step(y, t, dt);
                t += dt;
                rge.add_point(t, y);
            }
        }

        StabilityResult res = classify_stability(Mh, Mt);
        if (res.status != 2 && res.status != 3) {
            std::cerr << "Point is not metastable/unstable; no action curve." << std::endl;
            return 1;
        }
        double mu1 = res.mu_inst;

        // Scale of the minimal lambda_eff, for the search bracket center.
        double t_min = 2.0 * std::log(Mt);
        {
            double best_lam = 1.0;
            for (double tt = 2.0 * std::log(v); tt < 2.0 * std::log(planck_mass); tt += 0.1) {
                double lam = get_lambda_eff(rge.get_params(tt));
                if (lam < best_lam) { best_lam = lam; t_min = tt; }
            }
        }
        double R_opt = -1.0;
        find_minimum_action(rge, mu1, t_min, &R_opt);

        std::ofstream out(out_dir + "/action_curve.csv");
        if (!out.is_open()) { std::cerr << "cannot open action_curve.csv" << std::endl; return 1; }
        out << std::setprecision(12);
        out << "R_over_mu1,S_total,S_kinetic,S_potential\n";
        double logR_opt = std::log(R_opt);
        for (int i = 0; i <= 400; ++i) {
            double logR = logR_opt - 2.0 + 4.0 * i / 400.0;
            double R = std::exp(logR);
            ActionEvaluation e = evaluate_action_components(rge, mu1, R);
            // Skip the fenced region (lambda_R >= 0 or the kinetic shortcut):
            // genuine actions at these points are orders of magnitude below.
            if (e.total > 1e6) continue;
            out << R << "," << e.total << "," << e.kinetic << "," << e.potential << "\n";
        }
    }

    std::cout << "Diagnostics written to " << out_dir << std::endl;
    return 0;
}
