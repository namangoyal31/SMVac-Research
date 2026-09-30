#pragma once
#include <SMVacuumDecay/RGE.hpp>

namespace SMVacuumDecay {

// Production numerical settings (exposed so validation tests can vary them;
// see docs/validation.md for the measured convergence).
constexpr int kDefaultQuadraturePoints = 2048; // Simpson subintervals on u
constexpr double kDefaultRgeStep = 0.1;        // RK4 step in t = ln(mu^2)

// Decomposition of the ansatz action at fixed profile scale R.
struct ActionEvaluation {
    double total;     // S = S_kinetic + S_potential
    double kinetic;   // closed-form conformal kinetic term 16 pi^2 / (3 |lambda_R|)
    double potential; // Simpson integral of 2 pi^2 r^3 V(phi_R(r))
};

// RG-improved quartic coupling lambda_eff evaluated at the physical scale
// h = phi_dimless * mu_inst, capped at the Planck mass and frozen below the
// electroweak scale (where the high-field potential is not applicable).
double get_pure_sm_lambda(RGEHelper& rge, double phi_dimless, double mu_inst);

// Potential part of the ansatz action, integrand in the compactified
// variable u (see docs/numerical-method.md).
double integrand_u(RGEHelper& rge, double u, double R, double mu_inst, double prefactor);

// Action of the Fubini-Lipatov profile of scale R (in units of 1/mu_inst).
// quad_points overrides the production quadrature setting (for convergence
// tests).
ActionEvaluation evaluate_action_components(RGEHelper& rge, double mu_inst, double R,
                                            int quad_points = kDefaultQuadraturePoints);
double evaluate_action_at_R(RGEHelper& rge, double mu_inst, double R,
                            int quad_points = kDefaultQuadraturePoints);

// Minimum of the action over R, found by golden-section search on
// log(R) in a bracket around the naive conformal estimate. If R_opt_out is
// non-null it receives the optimal radius (in units of 1/mu_inst).
double find_minimum_action(RGEHelper& rge, double mu_inst, double t_min_lambda,
                           double* R_opt_out = nullptr);

// Full single-point stability analysis: status code, optimized actions,
// threshold and diagnostic scales.
struct StabilityResult {
    int status;        // 1 stable, 2 metastable, 3 unstable, 4 perturbativity lost
    double S_exact;    // optimized RG-improved ansatz action (-1 if stable)
    double S_approx;   // conformal estimate 8 pi^2 / (3 |lambda_min|) (-1 if stable)
    double S_kinetic;  // kinetic part of S_exact (-1 if stable)
    double S_potential;// potential part of S_exact (-1 if stable)
    double S_threshold;// metastability threshold 4 ln(mu1 * t_universe)
    double mu_inst;    // scale at which lambda_eff first crosses zero (-1 if stable)
    double lambda_min; // minimum of lambda_eff over the running range
    double Mh;         // Higgs mass input (echoed for convenience)
};

StabilityResult classify_stability(double Mh, double Mt,
                                   double rge_step = kDefaultRgeStep);

}
