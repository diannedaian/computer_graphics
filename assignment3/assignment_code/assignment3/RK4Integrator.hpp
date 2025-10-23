#ifndef RK4_INTEGRATOR_H_
#define RK4_INTEGRATOR_H_

#include "IntegratorBase.hpp"

namespace GLOO {
template <class TSystem, class TState>
class RK4Integrator : public IntegratorBase<TSystem, TState> {
  TState Integrate(const TSystem& system,
                   const TState& state,
                   float start_time,
                   float dt) const override {
    // RK4 integration: k1, k2, k3, k4
    auto k1 = system.ComputeTimeDerivative(state, start_time);
    auto k2 = system.ComputeTimeDerivative(state + dt * k1 * 0.5f, start_time + dt * 0.5f);
    auto k3 = system.ComputeTimeDerivative(state + dt * k2 * 0.5f, start_time + dt * 0.5f);
    auto k4 = system.ComputeTimeDerivative(state + dt * k3, start_time + dt);

    return state + dt * (k1 + 2.0f * k2 + 2.0f * k3 + k4) / 6.0f;
  }
};
}  // namespace GLOO

#endif
