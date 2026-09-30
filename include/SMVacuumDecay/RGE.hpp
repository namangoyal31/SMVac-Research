#pragma once

// Standard Model renormalization-group evolution for electroweak vacuum
// stability studies.
//
// Conventions (see docs/numerical-method.md for the full statement):
//   * Natural units, hbar = c = 1; all dimensional quantities in GeV.
//   * The Higgs field is parametrized by the real background phi = <h^0> and
//     the potential is written V = lambda*phi^4/4, so that at tree level
//     M_h^2 = 2*lambda*v^2.
//   * The hypercharge coupling g1 is GUT-normalized: g1 = sqrt(5/3)*g_Y.
//   * The RG evolution variable is t = ln(mu^2), and every beta function in
//     this library is dX/dt = dX/d ln(mu^2) = (1/2) dX/d ln(mu). The rates are
//     therefore HALF the standard d/d ln(mu) coefficients quoted in the
//     literature; the factor is absorbed consistently in rk4_single_step.
//
// The beta functions implement the 3-loop SM running as compiled in
// Buttazzo et al., JHEP 12 (2013) 089 [arXiv:1307.3536], with an additional
// 4-loop pure-gauge g3 contribution (see docs/audit-notes.md, finding 13).

#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>

namespace SMVacuumDecay {

// Physical constants (GeV unless stated otherwise).
const double pi = 3.14159265358979323846;
const double PI2 = pi * pi;
const double LOOP1 = 16.0 * PI2;    // (16 pi^2)     -- one-loop suppression
const double LOOP2 = LOOP1 * LOOP1; // (16 pi^2)^2   -- two-loop suppression
const double LOOP3 = LOOP2 * LOOP1; // (16 pi^2)^3   -- three-loop suppression
const double LOOP4 = LOOP2 * LOOP2; // (16 pi^2)^4   -- four-loop suppression

const double v = 246.22;            // Higgs vev, (sqrt(2) G_F)^(-1/2)
const double Mtau = 1.777;          // tau pole mass
const double Mb = 4.0;              // bottom mass used for yb(Mt)
const double alpha3_at_Mz = 0.1184; // alpha_s(M_Z), central value
const double planck_mass = 1.22e19; // non-reduced Planck mass

// Cosmology constants for the metastability criterion.
// The age of the universe is taken as a conservative 10 Gyr lower bound;
// the current 13.8 Gyr estimate would raise the action threshold
// S_threshold = 4 ln(mu * t_U) by only ~1.3 (see docs/theory.md).
const double universe_age_Gyr = 10.0;
const double seconds_per_year = 3.1536e7;  // 365 d
const double hbar_GeV_second = 6.582e-25;  // hbar [GeV s]
const double universe_age_GeV_minus1 =
    universe_age_Gyr * 1e9 * seconds_per_year / hbar_GeV_second;

struct StandardModelParameters {
    double g1, g2, g3, yt, yb, ytau, lambda, phi = 0;

    StandardModelParameters operator+(const StandardModelParameters& other) const {
        return {g1+other.g1, g2+other.g2, g3+other.g3, yt+other.yt, yb+other.yb, ytau+other.ytau, lambda+other.lambda, phi+other.phi};
    }
    StandardModelParameters operator*(double scalar) const {
        return {g1*scalar, g2*scalar, g3*scalar, yt*scalar, yb*scalar, ytau*scalar, lambda*scalar, phi*scalar};
    }
};

// NNLO electroweak matching conditions at mu = Mt, linearized around the
// central masses (Mt = 173.34 GeV), with the Mh-dependence carried by the
// tree-level relation lambda_tree = Mh^2/(2 v^2) plus a constant NNLO offset.
// The g1, g2 central values follow Buttazzo et al. (2013), eqs. (58)-(59).
// See docs/audit-notes.md (finding 16) for the accuracy of the linearization.
StandardModelParameters get_nnlo_matching(double Mh, double Mt);

// Beta functions in the convention dX/d ln(mu^2) described above.
double betaG1sq(const StandardModelParameters& p);
double betaG2sq(const StandardModelParameters& p);
double betaG3sq(const StandardModelParameters& p);
double betaLambda(const StandardModelParameters& p);
double betaYt2(const StandardModelParameters& p);
double betaYb2(const StandardModelParameters& p);
double betaYtau2(const StandardModelParameters& p);

// One fixed-step RK4 step of the SM RGE system over the evolution variable
// t = ln(mu^2).
StandardModelParameters rk4_single_step(const StandardModelParameters& y, double t, double dt);

// Tabulated RG trajectory with piecewise-linear interpolation in t.
// The table is filled by the caller at a fixed step size; the interpolation
// error is controlled by that step size (see docs/validation.md for the
// measured convergence).
class RGEHelper {
private:
    std::vector<std::pair<double, StandardModelParameters>> table;
    bool sorted = false;
public:
    void clear();
    void add_point(double t, const StandardModelParameters& p);
    void sort_table();
    // Linear interpolation in t; outside the tabulated range the parameters
    // are frozen at the nearest endpoint.
    StandardModelParameters get_params(double t);
};

} // namespace SMVacuumDecay
