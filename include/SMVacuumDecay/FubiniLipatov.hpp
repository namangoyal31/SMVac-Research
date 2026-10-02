#pragma once
#include <tuple>

namespace SMVacuumDecay {

// Strict conformal (Fubini-Lipatov) estimate evaluated on the RG-evolved
// couplings: the RG equations are run up to the Planck scale, the minimum
// lambda_eff is located, and S_conformal = 8 pi^2 / (3 |lambda_min|) is
// formed from it. Returns (status, S_conformal, Mh) with the same status
// codes as classify_stability (1 stable, 2 metastable, 3 unstable,
// 4 perturbativity lost). This is the standard leading-order stability
// estimate used in the literature (e.g. Isidori-Ridolfi-Strumia); it is not
// a bounce computation and does not reproduce any part of Buttazzo et al.'s
// NNLO analysis despite the historical function name.
std::tuple<int, double, double> classify_conformal(double Mh_input, double Mt);

}
