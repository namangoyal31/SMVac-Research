// Regenerates the frozen reference tables used by the regression tests:
//   reference/v1.1/rge_reference.csv          -- RG trajectory for the SM point
//   reference/v1.1/potential_reference.csv    -- potential along the trajectory
//   reference/v1.1/benchmark_results*.json    -- single-point benchmark values
//
// The tables describe the benchmark point Mh = 125.1 GeV, Mt = 173.1 GeV and
// the low-Higgs point Mh = 5 GeV, Mt = 105 GeV, and are produced with the
// same fixed-step integrator settings (dt = 0.1 in t = ln(mu^2)) used by the
// production pipeline. Re-run this tool whenever the physics or the
// numerical settings change intentionally, and record the change in
// docs/audit-notes.md and CHANGELOG.md.

#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <SMVacuumDecay/FubiniLipatov.hpp>
#include <SMVacuumDecay/EffectivePotential.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#endif

using namespace SMVacuumDecay;

namespace {

void ensure_dir(const std::string& path) {
#ifdef _WIN32
    _mkdir(path.c_str());
#else
    mkdir(path.c_str(), 0755);
#endif
}

void write_rge_table(const std::string& filename, double Mh, double Mt) {
    std::ofstream out(filename);
    if (!out.is_open()) { std::cerr << "cannot open " << filename << "\n"; std::exit(1); }
    out << "scale_mu,g1,g2,g3,yt,lambda\n";

    StandardModelParameters y = get_nnlo_matching(Mh, Mt);
    double t_match = 2 * std::log(Mt);
    double t_start = 2 * std::log(1.0);
    double t_planck = 2 * std::log(planck_mass);
    double dt = 0.1;

    // Backward from the matching scale down to mu = 1 GeV, then forward up to
    // the Planck scale; identical to the trajectory built by classify_stability
    // plus the low-scale extension.
    StandardModelParameters y_rev = y;
    double t_rev = t_match;
    std::vector<std::pair<double, StandardModelParameters>> rows;
    rows.push_back({t_rev, y_rev});
    while (t_rev > t_start) {
        if (t_rev - dt < t_start) dt = t_rev - t_start;
        y_rev = rk4_single_step(y_rev, t_rev, -dt);
        t_rev -= dt;
        rows.push_back({t_rev, y_rev});
    }
    StandardModelParameters y_fwd = y;
    double t_fwd = t_match;
    dt = 0.1;
    while (t_fwd < t_planck) {
        if (t_fwd + dt > t_planck) dt = t_planck - t_fwd;
        y_fwd = rk4_single_step(y_fwd, t_fwd, dt);
        t_fwd += dt;
        rows.push_back({t_fwd, y_fwd});
    }

    out << std::setprecision(17);
    for (const auto& row : rows) {
        double mu = std::exp(row.first / 2.0);
        const auto& p = row.second;
        out << mu << "," << p.g1 << "," << p.g2 << "," << p.g3 << ","
            << p.yt << "," << p.lambda << "\n";
    }
}

void write_potential_table(const std::string& filename, double Mh, double Mt) {
    std::ofstream out(filename);
    if (!out.is_open()) { std::cerr << "cannot open " << filename << "\n"; std::exit(1); }
    out << "phi,lambda_eff,V\n";

    StandardModelParameters y = get_nnlo_matching(Mh, Mt);
    double t_match = 2 * std::log(Mt);
    double t_start = 2 * std::log(1.0);
    double t_planck = 2 * std::log(planck_mass);
    double dt = 0.1;

    RGEHelper rge;
    rge.add_point(t_match, y);
    StandardModelParameters y_rev = y;
    double t_rev = t_match;
    while (t_rev > t_start) {
        if (t_rev - dt < t_start) dt = t_rev - t_start;
        y_rev = rk4_single_step(y_rev, t_rev, -dt);
        t_rev -= dt;
        rge.add_point(t_rev, y_rev);
    }
    StandardModelParameters y_fwd = y;
    double t_fwd = t_match;
    dt = 0.1;
    while (t_fwd < t_planck) {
        if (t_fwd + dt > t_planck) dt = t_planck - t_fwd;
        y_fwd = rk4_single_step(y_fwd, t_fwd, dt);
        t_fwd += dt;
        rge.add_point(t_fwd, y_fwd);
    }

    // phi = v * exp(0.1 k), k = 0, 1, ... up to the Planck scale.
    out << std::setprecision(17);
    double phi = v;
    while (phi < planck_mass) {
        double t = 2.0 * std::log(phi);
        StandardModelParameters p = rge.get_params(t);
        double lam_eff = get_lambda_eff(p);
        out << phi << "," << lam_eff << "," << V_eff(phi, rge) << "\n";
        phi = phi * std::exp(0.1);
    }
}

struct BenchmarkPoint { const char* name; double Mh, Mt; };

void write_benchmark_json(const std::string& filename,
                          const std::vector<BenchmarkPoint>& points) {
    std::ofstream out(filename);
    if (!out.is_open()) { std::cerr << "cannot open " << filename << "\n"; std::exit(1); }
    out << std::setprecision(17) << "{";
    bool first = true;
    for (const auto& pt : points) {
        auto fl = classify_conformal(pt.Mh, pt.Mt);
        auto num = classify_stability(pt.Mh, pt.Mt);
        if (!first) out << ", ";
        first = false;
        out << "\"" << pt.name << "\": {\"Mh\": " << pt.Mh << ", \"Mt\": " << pt.Mt
            << ", \"Conformal_action\": " << std::get<1>(fl)
            << ", \"Trial_action\": " << num.S_trial
            << ", \"Status\": " << num.status
            << ", \"Method_flag\": " << num.method_flag << "}";
    }
    out << "}\n";
}

void write_benchmark_json_analytical(const std::string& filename,
                                     const std::vector<BenchmarkPoint>& points) {
    std::ofstream out(filename);
    if (!out.is_open()) { std::cerr << "cannot open " << filename << "\n"; std::exit(1); }
    out << std::setprecision(17) << "{";
    bool first = true;
    for (const auto& pt : points) {
        auto fl = classify_conformal(pt.Mh, pt.Mt);
        if (!first) out << ", ";
        first = false;
        out << "\"" << pt.name << "\": {\"Conformal_action\": " << std::get<1>(fl)
            << ", \"Status\": " << std::get<0>(fl) << "}";
    }
    out << "}\n";
}

void write_benchmark_json_numerical(const std::string& filename,
                                    const std::vector<BenchmarkPoint>& points) {
    std::ofstream out(filename);
    if (!out.is_open()) { std::cerr << "cannot open " << filename << "\n"; std::exit(1); }
    out << std::setprecision(17) << "{";
    bool first = true;
    for (const auto& pt : points) {
        auto num = classify_stability(pt.Mh, pt.Mt);
        if (!first) out << ", ";
        first = false;
        out << "\"" << pt.name << "\": {\"Trial_action\": " << num.S_trial
            << ", \"Status\": " << num.status
            << ", \"Method_flag\": " << num.method_flag << "}";
    }
    out << "}\n";
}

} // namespace

int main() {
    ensure_dir("reference");
    ensure_dir("reference/v1.1");

    write_rge_table("reference/v1.1/rge_reference.csv", 125.1, 173.1);
    // The potential table uses the low-Higgs point, matching the original
    // reference set (a second, deeply metastable region of parameter space).
    write_potential_table("reference/v1.1/potential_reference.csv", 5.0, 105.0);

    std::vector<BenchmarkPoint> points = {
        {"SM125", 125.1, 173.1},
        {"LowMh", 5.0, 105.0},
    };
    write_benchmark_json("reference/v1.1/benchmark_results.json", points);
    write_benchmark_json_analytical("reference/v1.1/benchmark_results_analytical.json", points);
    write_benchmark_json_numerical("reference/v1.1/benchmark_results_numerical.json", points);

    std::cout << "Reference tables written to reference/v1.1/" << std::endl;
    return 0;
}
