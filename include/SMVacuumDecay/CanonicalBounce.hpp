#pragma once
#include <SMVacuumDecay/RGE.hpp>

namespace SMVacuumDecay {

// Production numerical settings (exposed so validation tests can vary them;
// see docs/validation.md for the measured convergence).
constexpr int kDefaultQuadraturePoints = 2048; // Simpson subintervals on u
constexpr double kDefaultRgeStep = 0.1;        // RK4 step in t = ln(mu^2)
constexpr double kDefaultC6 = 1.0;             // Planck-suppressed phi^6 coefficient

// Distinguishes how the trial action was obtained at a given point, so that
// method failures are never read as physical classifications:
enum MethodFlag : int {
    MethodOK = 0,                // trial action evaluated normally
    MethodKineticShortcut = 1,   // action = closed-form kinetic term only
                                 // (|lambda_R| <~ 3.3e-4): rough overestimate
    MethodFence = 2,             // no radius in the search bracket gives
                                 // lambda_R < 0: NO action computed (S_trial = NaN)
    MethodAnsatzFailed = 3,      // trial action evaluated but non-positive:
                                 // ansatz breakdown, classification undefined
    MethodPerturbativityLost = 4,// RGE left perturbative range before M_Pl
};

// Physical status of a grid point. status = 0 means "undetermined" -- the
// method_flag then records why (Fence, AnsatzFailed); those points must not
// be counted as metastable or unstable.
//   1 stable (lambda_eff >= 0 up to M_Pl), 2 metastable, 3 unstable,
//   4 perturbativity lost, 0 undetermined.
struct StabilityResult {
    int status;        // 0/1/2/3/4 as described above
    int method_flag;   // MethodFlag value
    double S_trial;    // optimized trial-profile action (NaN if not computed;
                       // may be negative where the ansatz breaks down)
    double S_conformal;// strict conformal estimate 8 pi^2 / (3 |lambda_min|)
    double S_kinetic;  // kinetic part of S_trial at the optimal radius
    double S_potential;// potential part of S_trial at the optimal radius
    double S_threshold;// metastability threshold 4 ln(mu1 * t_universe)
    double mu_inst;    // scale at which lambda_eff first crosses zero (-1 if stable)
    double lambda_min; // minimum of lambda_eff over the running range
    double Mh;         // Higgs mass input (echoed for convenience)
};

// Decomposition of the trial action at fixed profile scale R.
struct ActionEvaluation {
    double total;     // S = S_kinetic + S_potential
    double kinetic;   // closed-form conformal kinetic term 16 pi^2 / (3 |lambda_R|)
    double potential; // Simpson integral of 2 pi^2 r^3 V(phi_R(r))
};

// RG-improved quartic coupling lambda_eff evaluated at the physical scale
// h = phi_dimless * mu_inst, capped at the Planck mass and frozen below the
// electroweak scale (where the high-field potential is not applicable).
double get_pure_sm_lambda(RGEHelper& rge, double phi_dimless, double mu_inst);

// Potential part of the trial action, integrand in the compactified
// variable u (see docs/numerical-method.md). c6 is the coefficient of the
// Planck-suppressed phi^6 operator (docs/limitations.md, item 3).
double integrand_u(RGEHelper& rge, double u, double R, double mu_inst,
                   double prefactor, double c6);

// Trial action of the Fubini-Lipatov profile of scale R (in units of
// 1/mu_inst). quad_points overrides the production quadrature setting and c6
// the Planck-operator coefficient (both for sensitivity studies).
ActionEvaluation evaluate_action_components(RGEHelper& rge, double mu_inst, double R,
                                            int quad_points = kDefaultQuadraturePoints,
                                            double c6 = kDefaultC6);
double evaluate_action_at_R(RGEHelper& rge, double mu_inst, double R,
                            int quad_points = kDefaultQuadraturePoints,
                            double c6 = kDefaultC6);

// Minimum of the trial action over R, found by golden-section search on
// log(R) in a bracket around the naive conformal estimate (bracket spans
// +-8 in ln R, i.e. ~7 decades in R). If R_opt_out is non-null it receives
// the optimal radius (in units of 1/mu_inst). Returns >= 1e100 if no radius
// in the bracket gives lambda_R < 0.
double find_minimum_action(RGEHelper& rge, double mu_inst, double t_min_lambda,
                           double* R_opt_out = nullptr, double c6 = kDefaultC6);

// Full single-point analysis: physical status, method flag, optimized
// trial-profile actions, threshold and diagnostic scales. c6 is the
// coefficient of the Planck-suppressed phi^6 operator.
StabilityResult classify_stability(double Mh, double Mt,
                                   double rge_step = kDefaultRgeStep,
                                   double c6 = kDefaultC6);

}
