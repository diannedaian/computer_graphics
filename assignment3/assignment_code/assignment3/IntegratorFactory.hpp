#ifndef INTEGRATOR_FACTORY_H_
#define INTEGRATOR_FACTORY_H_

#include "IntegratorBase.hpp"
#include "ForwardEulerIntegrator.hpp"
#include "TrapezoidalIntegrator.hpp"
#include "RK4Integrator.hpp"

#include <stdexcept>

#include "gloo/utils.hpp"

#include "IntegratorType.hpp"

namespace GLOO {
class IntegratorFactory {
 public:
  template <class TSystem, class TState>
  static std::unique_ptr<IntegratorBase<TSystem, TState>> CreateIntegrator(
      IntegratorType type) {
    switch (type) {
      case IntegratorType::Euler:
        return make_unique<ForwardEulerIntegrator<TSystem, TState>>();
      case IntegratorType::Trapezoidal:
        return make_unique<TrapezoidalIntegrator<TSystem, TState>>();
      case IntegratorType::RK4:
        return make_unique<RK4Integrator<TSystem, TState>>();
      default:
        throw std::runtime_error("Unknown integrator type");
    }
  }
};
}  // namespace GLOO

#endif
