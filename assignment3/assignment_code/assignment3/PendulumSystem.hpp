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
  float drag_k = 1.0f;  // Increased drag for better stability with cloth
  glm::vec3 gravity = glm::vec3(0.0f, -9.8f, 0.0f);

  // Wind properties
  bool wind_enabled = false;
  glm::vec3 wind_direction = glm::vec3(1.0f, 0.0f, 0.0f);  // Default wind direction
  float wind_strength = 5.0f;  // Wind force magnitude
  float wind_frequency = 0.5f;  // How often wind changes
  float wind_time = 0.0f;  // Time accumulator for wind variation

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

  // ---- Wind control methods ----
  void ToggleWind() { wind_enabled = !wind_enabled; }
  void SetWindEnabled(bool enabled) { wind_enabled = enabled; }
  bool IsWindEnabled() const { return wind_enabled; }

  // Update wind direction with time-based variation
  void UpdateWind(float dt) {
    if (!wind_enabled) return;

    wind_time += dt;

    // Create gentle, time-varying wind
    float base_angle = wind_time * wind_frequency;
    float variation = 0.3f * sin(wind_time * 0.7f) + 0.2f * sin(wind_time * 1.3f);

    wind_direction.x = cos(base_angle + variation);
    wind_direction.y = 0.1f * sin(base_angle + variation);  // Slight vertical component
    wind_direction.z = 0.2f * sin(wind_time * 0.5f);  // Some z-direction variation

    // Normalize the direction
    wind_direction = glm::normalize(wind_direction);
  }

  // Get current wind force for a particle
  glm::vec3 GetWindForce(int particle_index) const {
    if (!wind_enabled || fixed[particle_index]) {
      return glm::vec3(0.0f);
    }

    // Add some randomness based on particle position for more realistic effect
    float noise = 0.5f + 0.5f * sin(particle_index * 0.1f + wind_time * 2.0f);
    return wind_strength * wind_direction * noise;
  }

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

      // Add wind force
      force += GetWindForce(i);

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
