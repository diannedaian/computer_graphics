#pragma once
#include "gloo/SceneNode.hpp"
#include "IntegratorFactory.hpp"
#include "PendulumSystem.hpp"

#include "gloo/shaders/SimpleShader.hpp"
#include "gloo/debug/PrimitiveFactory.hpp"
#include "gloo/components/ShadingComponent.hpp"
#include "gloo/components/RenderingComponent.hpp"
#include "gloo/InputManager.hpp"

namespace GLOO {

class PendulumSystemNode : public SceneNode {
 public:
  PendulumSystemNode(IntegratorType integrator_type, float dt, int num_particles)
      : m_dt(dt), m_time(0.0f), m_num_particles(num_particles), m_accumulator(0.0) {
    std::cerr << "Initializing PendulumSystemNode..." << std::endl;

    // ---- 1.  Build the pendulum system  ----
    const float rest_len = 0.3f;  // Shorter segments for better stability
    const float k_spring = 50.0f;  // Reduced spring constant for stability

    for (int i = 0; i < m_num_particles; i++) {
      glm::vec3 pos(0.0f, -i * rest_len, 0.0f);
      m_system.AddParticle(pos, 1.0f, glm::vec3(0.0f));
    }

    for (int i = 0; i < m_num_particles - 1; i++) {
      m_system.AddSpring(i, i + 1, k_spring, rest_len);
    }

    // Fix the first particle (the top)
    m_system.FixParticle(0);

    // Keep a copy of the starting state
    m_state = m_system.GetInitialState();

    // ---- 2.  Create the integrator  ----
    m_integrator = IntegratorFactory::CreateIntegrator<PendulumSystem, ParticleState>(integrator_type);

    // ---- 3.  Build visual spheres ----
    auto shader = std::make_shared<SimpleShader>();
    auto sphere = PrimitiveFactory::CreateSphere(0.03f, 25, 25);
    auto sphere_shared = std::shared_ptr<VertexObject>(sphere.release());

    for (int i = 0; i < m_num_particles; i++) {
      auto child = make_unique<SceneNode>();
      child->CreateComponent<ShadingComponent>(shader);
      child->CreateComponent<RenderingComponent>(sphere_shared);
      child->GetTransform().SetPosition(m_state.positions[i]);
      SceneNode* ptr = child.get();
      AddChild(std::move(child));
      m_particle_nodes.push_back(ptr);
    }

    std::cerr << "PendulumSystemNode initialized with "
              << m_num_particles << " particles." << std::endl;
  }

  // ---- 4.  Per-frame update ----
  void Update(double delta_time) override {
    if (InputManager::GetInstance().IsKeyPressed('R'))
      Reset();

    // Accumulate time so simulation is frame-rate independent
    m_accumulator += delta_time;
    while (m_accumulator >= m_dt) {
      m_state = m_integrator->Integrate(m_system, m_state, m_time, m_dt);
      m_time += m_dt;
      m_accumulator -= m_dt;
    }

    // Update visuals
    for (int i = 0; i < m_num_particles; i++) {
      m_particle_nodes[i]->GetTransform().SetPosition(m_state.positions[i]);
    }
  }

  // ---- 5.  Reset key ----
  void Reset() {
    m_state = m_system.GetInitialState();
    m_time = 0.0f;
    m_accumulator = 0.0;
    for (int i = 0; i < m_num_particles; i++) {
      m_particle_nodes[i]->GetTransform().SetPosition(m_state.positions[i]);
    }
  }

 private:
  float m_dt, m_time;
  double m_accumulator;
  int m_num_particles;

  ParticleState m_state;
  PendulumSystem m_system;

  std::unique_ptr<IntegratorBase<PendulumSystem, ParticleState>> m_integrator;
  std::vector<SceneNode*> m_particle_nodes;
};

}  // namespace GLOO
