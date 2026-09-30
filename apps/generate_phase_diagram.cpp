// Parameter-scan application: classifies the (Mh, Mt) plane point by point
// with either the analytical (strict conformal) or the numerical
// (RG-improved ansatz) estimator, writing one CSV per chunk.
//
// Usage:
//   generate_phase_diagram --analytical|--numerical
//                          [--mt-min A] [--mt-max B] [--mh-min C] [--mh-max D]
//                          [--step s] [--start i] [--end j] [--output-dir DIR]
//
// The grid is enumerated Mt-major (Mt outer, Mh inner); [start, end) selects
// an index range for chunked runs. Default grid covers the physically
// relevant region Mt in [160, 185], Mh in [110, 140] GeV with 0.25 GeV steps.

#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <SMVacuumDecay/FubiniLipatov.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>

using namespace SMVacuumDecay;

namespace {

struct ScanConfig {
    bool use_analytical = false;
    double mt_min = 160.0, mt_max = 185.0;
    double mh_min = 110.0, mh_max = 140.0;
    double step = 0.25;
    size_t start_idx = 0;
    size_t end_idx = 0;      // 0 -> full range
    std::string output_dir = "results";
};

double parse_arg(const char* name, char* value) {
    char* end = nullptr;
    double v = std::strtod(value, &end);
    if (end == value) {
        std::cerr << "Invalid numeric value for " << name << ": " << value << std::endl;
        std::exit(1);
    }
    return v;
}

ScanConfig parse_args(int argc, char* argv[]) {
    ScanConfig cfg;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        auto next_value = [&]() -> char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << arg << std::endl;
                std::exit(1);
            }
            return argv[++i];
        };
        if (arg == "--analytical") cfg.use_analytical = true;
        else if (arg == "--numerical") cfg.use_analytical = false;
        else if (arg == "--mt-min") cfg.mt_min = parse_arg(arg.c_str(), next_value());
        else if (arg == "--mt-max") cfg.mt_max = parse_arg(arg.c_str(), next_value());
        else if (arg == "--mh-min") cfg.mh_min = parse_arg(arg.c_str(), next_value());
        else if (arg == "--mh-max") cfg.mh_max = parse_arg(arg.c_str(), next_value());
        else if (arg == "--step") cfg.step = parse_arg(arg.c_str(), next_value());
        else if (arg == "--start") cfg.start_idx = static_cast<size_t>(parse_arg(arg.c_str(), next_value()));
        else if (arg == "--end") cfg.end_idx = static_cast<size_t>(parse_arg(arg.c_str(), next_value()));
        else if (arg == "--output-dir") cfg.output_dir = next_value();
        else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            std::exit(1);
        }
    }
    if (cfg.step <= 0) { std::cerr << "--step must be positive" << std::endl; std::exit(1); }
    if (cfg.mt_max < cfg.mt_min || cfg.mh_max < cfg.mh_min) {
        std::cerr << "Empty scan range" << std::endl;
        std::exit(1);
    }
    return cfg;
}

} // namespace

int main(int argc, char* argv[]) {
    ScanConfig cfg = parse_args(argc, argv);

    // Mt-major enumeration: n_mt rows of n_mh points each.
    size_t n_mt = static_cast<size_t>(std::lround((cfg.mt_max - cfg.mt_min) / cfg.step)) + 1;
    size_t n_mh = static_cast<size_t>(std::lround((cfg.mh_max - cfg.mh_min) / cfg.step)) + 1;
    size_t total = n_mt * n_mh;

    size_t end_idx = cfg.end_idx;
    if (end_idx == 0 || end_idx > total) end_idx = total;
    if (cfg.start_idx >= end_idx) {
        std::cerr << "Empty index range [" << cfg.start_idx << ", " << end_idx << ")" << std::endl;
        return 1;
    }

    std::string prefix = cfg.use_analytical ? "analytical" : "numerical";
    std::string filename = cfg.output_dir + "/" + prefix + "_chunk_" + std::to_string(cfg.start_idx) + ".csv";

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "ERROR: could not open output file '" << filename
                  << "'. Does the output directory exist? The Python runners "
                  << "create it automatically; otherwise create it manually."
                  << std::endl;
        return 1;
    }

    if (cfg.use_analytical) {
        file << "Mt,Mh_calc,Stability,S_approx\n";
    } else {
        file << "Mt,Mh_calc,Stability,S_exact,S_approx,S_kinetic,S_potential,S_threshold,mu_inst,lambda_min\n";
    }

    std::cout << "Analyzing " << prefix << " points [" << cfg.start_idx << ", " << end_idx
              << ") of " << total << " (Mt in [" << cfg.mt_min << ", " << cfg.mt_max
              << "], Mh in [" << cfg.mh_min << ", " << cfg.mh_max << "], step "
              << cfg.step << ")" << std::endl;

    for (size_t i = cfg.start_idx; i < end_idx; ++i) {
        size_t i_mt = i / n_mh;
        size_t i_mh = i % n_mh;
        double Mt = cfg.mt_min + cfg.step * i_mt;
        double Mh_input = cfg.mh_min + cfg.step * i_mh;

        if (cfg.use_analytical) {
            auto res = classify_buttazzo(Mh_input, Mt);
            file << Mt << "," << std::get<2>(res) << "," << std::get<0>(res)
                 << "," << std::get<1>(res) << "\n";
        } else {
            StabilityResult res = classify_stability(Mh_input, Mt);
            file << Mt << "," << res.Mh << "," << res.status
                 << "," << res.S_exact << "," << res.S_approx
                 << "," << res.S_kinetic << "," << res.S_potential
                 << "," << res.S_threshold << "," << res.mu_inst
                 << "," << res.lambda_min << "\n";
        }

        if (i % 1000 == 0) std::cout << "Processed " << i << std::endl;
    }

    bool write_ok = file.good();
    file.close();
    if (!write_ok) {
        std::cerr << "ERROR: writing '" << filename << "' failed partway through." << std::endl;
        return 1;
    }
    std::cout << "Done. Saved to " << filename << std::endl;
    return 0;
}
