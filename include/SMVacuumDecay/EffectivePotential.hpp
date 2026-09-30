#pragma once
#include <SMVacuumDecay/RGE.hpp>

namespace SMVacuumDecay {

// RG-improved 1-loop effective quartic coupling:
//   lambda_eff(phi) = lambda(mu=phi) + Delta_lambda_CW(couplings at mu=phi),
// with the Coleman-Weinberg correction from the top quark and the W/Z gauge
// bosons in the high-field approximation (Goldstone and Higgs-loop terms,
// which carry no phi^4 weight, are omitted). Valid for phi >> v.
double get_lambda_eff(const StandardModelParameters& p);

// RG-improved effective potential in the high-field approximation,
//   V(phi) = lambda_eff(phi) * phi^4 / 4,
// with couplings interpolated from the tabulated RG trajectory at
// mu = phi. The tree-level mass term is deliberately absent: the conformal
// ansatz treats phi = 0 as the false vacuum (see docs/theory.md). Returns 0
// for phi <= 0.
double V_eff(double phi, RGEHelper& rge);

}
