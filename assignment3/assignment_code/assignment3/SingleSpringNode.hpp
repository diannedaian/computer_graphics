#pragma once
#include "gloo/SceneNode.hpp"
#include "PendulumSystem.hpp"
#include "IntegratorFactory.hpp"
#include "gloo/shaders/SimpleShader.hpp"
#include "gloo/debug/PrimitiveFactory.hpp"
#include "gloo/components/ShadingComponent.hpp"
#include "gloo/components/RenderingComponent.hpp"
#include "gloo/InputManager.hpp"

namespace GLOO {

class SingleSpringNode : public SceneNode {
 public:
  SingleSpringNode(IntegratorType integrator_type, float dt)
      : m_dt(dt), m_time(0.0f), m_accumulator(0.0) {
    std::cerr << "Initializing single-particle spring test..." << std::endl;

    // Create PendulumSystem with 2 particles: one fixed anchor, one moving mass
    m_system.AddParticle(glm::vec3(0.0f, 0.0f, 0.0f), 1.0f);   // anchor (fixed)
    m_system.AddParticle(glm::vec3(0.0f, -0.5f, 0.0f), 1.0f);  // moving particle
    m_system.FixParticle(0);

    // Connect them with a spring
    m_system.AddSpring(0, 1, 100.0f, 0.5f);
    m_state = m_system.GetInitialState();

    // Choose integrator
    m_integrator =
        IntegratorFactory::CreateIntegrator<PendulumSystem, ParticleState>(
            integrator_type);

    // Visuals: spheres for anchor and particle
    auto shader = std::make_shared<SimpleShader>();
    auto sphere = PrimitiveFactory::CreateSphere(0.03f, 25, 25);
    auto sphere_shared = std::shared_ptr<VertexObject>(sphere.release());

    for (int i = 0; i < 2; i++) {
      auto child = make_unique<SceneNode>();
      child->CreateComponent<ShadingComponent>(shader);
      child->CreateComponent<RenderingComponent>(sphere_shared);
      child->GetTransform().SetPosition(m_state.positions[i]);
      SceneNode* ptr = child.get();
      AddChild(std::move(child));
      m_particle_nodes.push_back(ptr);
    }

    std::cerr << "SingleSpringNode initialized." << std::endl;
  }

  void Update(double delta_time) override {
    if (InputManager::GetInstance().IsKeyPressed('R')) Reset();

    m_accumulator += delta_time;
    while (m_accumulator >= m_dt) {
      m_state = m_integrator->Integrate(m_system, m_state, m_time, m_dt);
      m_time += m_dt;
      m_accumulator -= m_dt;
    }

    for (int i = 0; i < 2; i++)
      m_particle_nodes[i]->GetTransform().SetPosition(m_state.positions[i]);
  }

  void Reset() {
    m_state = m_system.GetInitialState();
    m_time = 0.0f;
    m_accumulator = 0.0;
    for (int i = 0; i < 2; i++)
      m_particle_nodes[i]->GetTransform().SetPosition(m_state.positions[i]);
  }

 private:
  float m_dt, m_time;
  double m_accumulator;
  ParticleState m_state;
  PendulumSystem m_system;
  std::unique_ptr<IntegratorBase<PendulumSystem, ParticleState>> m_integrator;
  std::vector<SceneNode*> m_particle_nodes;
};

}  // namespace GLOO
