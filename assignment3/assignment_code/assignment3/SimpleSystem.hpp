#pragma once
#include "ParticleSystemBase.hpp"
#include <glm/glm.hpp>

namespace GLOO {

class SimpleSystem : public ParticleSystemBase {
public:
    ParticleState ComputeTimeDerivative(const ParticleState& state, float time) const override {
        ParticleState derivative;
        derivative.positions.resize(1);
        derivative.velocities.resize(1);

        glm::vec3 p = state.positions[0];
        derivative.positions[0] = glm::vec3(-p.y, p.x, 0.0f);
        derivative.velocities[0] = glm::vec3(0.0f);

        return derivative;
    }
};

}  // namespace GLOO
