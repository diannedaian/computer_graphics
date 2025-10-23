#pragma once
#include "ParticleSystemBase.hpp"
#include "ParticleState.hpp"
#include <glm/glm.hpp>
#include <vector>

namespace GLOO {

struct Spring {
  int i, j;          // particle indices
  float rest_length;
  float k;           // stiffness
};

class PendulumSystem : public ParticleSystemBase {
public:
  float drag_k = 0.5f;  // Increased drag for better stability
  glm::vec3 gravity = glm::vec3(0.0f, -9.8f, 0.0f);

  // Data storage
  std::vector<float> masses;
  std::vector<bool> fixed;
  std::vector<Spring> springs;
  ParticleState initial_state_;

  // ---- Add a new particle ----
  int AddParticle(const glm::vec3& pos, float mass,
                  const glm::vec3& vel = glm::vec3(0.0f)) {
    int idx = masses.size();
    masses.push_back(mass);
    fixed.push_back(false);
    initial_state_.positions.push_back(pos);
    initial_state_.velocities.push_back(vel);
    return idx;
  }

  // ---- Add spring between two particles ----
  void AddSpring(int i, int j, float k, float rest_len) {
    springs.push_back({i, j, rest_len, k});
  }

  // ---- Fix a particle (pin it) ----
  void FixParticle(int i) { fixed[i] = true; }

  // ---- Retrieve initial state (for Reset / start-up) ----
  const ParticleState& GetInitialState() const { return initial_state_; }

  // ---- Compute time derivative (core physics) ----
  ParticleState ComputeTimeDerivative(const ParticleState& state, float time) const override {
    ParticleState deriv;
    int n = state.positions.size();
    deriv.positions.resize(n);
    deriv.velocities.resize(n);

    for (int i = 0; i < n; i++) {
      deriv.positions[i] = state.velocities[i];
      deriv.velocities[i] = glm::vec3(0.0f);
    }

    // Compute forces for each particle
    for (int i = 0; i < n; i++) {
      if (fixed[i]) continue;

      glm::vec3 force = masses[i] * gravity;
      force += -drag_k * state.velocities[i]; // viscous drag

      // Sum spring forces
      for (const auto& s : springs) {
        if (s.i == i || s.j == i) {
          int j = (s.i == i) ? s.j : s.i;
          glm::vec3 d = state.positions[i] - state.positions[j];
          float dist = glm::length(d);
          if (dist > 1e-6f) {
            glm::vec3 dir = d / dist;
            force += -s.k * (dist - s.rest_length) * dir;
          }
        }
      }

      deriv.velocities[i] = force / masses[i];
    }

    // Keep fixed particles stationary
    for (int i = 0; i < n; i++) {
      if (fixed[i]) {
        deriv.positions[i] = glm::vec3(0.0f);
        deriv.velocities[i] = glm::vec3(0.0f);
      }
    }

    return deriv;
  }
};

}  // namespace GLOO
