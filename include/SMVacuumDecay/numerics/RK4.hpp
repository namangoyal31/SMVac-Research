#pragma once
#include <cmath>

namespace SMVacuumDecay {
namespace numerics {

// One classical 4th-order Runge-Kutta step for an arbitrary state type that
// provides operator+ (State, State) and operator* (State, double).
template <typename State, typename DerivFunc>
State rk4_single_step(const State& y, double t, double dt, DerivFunc dydt) {
    State k1 = dydt(y, t);
    State k2 = dydt(y + k1 * (0.5 * dt), t + 0.5 * dt);
    State k3 = dydt(y + k2 * (0.5 * dt), t + 0.5 * dt);
    State k4 = dydt(y + k3 * dt, t + dt);

    return y + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
}

} // namespace numerics
} // namespace SMVacuumDecay
