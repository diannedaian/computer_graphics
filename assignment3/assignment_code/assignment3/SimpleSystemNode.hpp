#pragma once
#include "gloo/SceneNode.hpp"
#include "IntegratorFactory.hpp"
#include "SimpleSystem.hpp"

#include "gloo/shaders/SimpleShader.hpp"
#include "gloo/debug/PrimitiveFactory.hpp"
#include "gloo/components/ShadingComponent.hpp"
#include "gloo/components/RenderingComponent.hpp"
#include "gloo/InputManager.hpp"

namespace GLOO {

class SimpleSystemNode : public SceneNode {
 public:
  SimpleSystemNode(IntegratorType integrator_type, float dt)
      : m_dt(dt), m_time(0.0f), m_accumulator(0.0) {
    std::cerr << "Initializing SimpleSystemNode..." << std::endl;

    // Initialize the particle state with a starting position
    m_state.positions.push_back(glm::vec3(1.0f, 0.0f, 0.0f));  // Start at (1, 0, 0)
    m_state.velocities.push_back(glm::vec3(0.0f, 0.0f, 0.0f)); // No initial velocity

    // Create the integrator
    m_integrator = IntegratorFactory::CreateIntegrator<SimpleSystem, ParticleState>(integrator_type);

    // Create visual sphere for the particle
    auto shader = std::make_shared<SimpleShader>();
    auto sphere = PrimitiveFactory::CreateSphere(0.05f, 25, 25);
    auto sphere_shared = std::shared_ptr<VertexObject>(sphere.release());

    auto particle_node = make_unique<SceneNode>();
    particle_node->CreateComponent<ShadingComponent>(shader);
    particle_node->CreateComponent<RenderingComponent>(sphere_shared);
    particle_node->GetTransform().SetPosition(m_state.positions[0]);
    m_particle_node = particle_node.get();
    AddChild(std::move(particle_node));

    // Create a trail to show the particle's path
    CreateTrail();

    std::cerr << "SimpleSystemNode initialized." << std::endl;
  }

  // Create a trail to visualize the particle's path
  void CreateTrail() {
    auto line_shader = std::make_shared<SimpleShader>();

    // Initialize trail positions
    m_trail_positions.clear();
    m_trail_positions.push_back(m_state.positions[0]);

    // Create trail line
    auto trail = make_unique<VertexObject>();
    auto positions = make_unique<PositionArray>(m_trail_positions);
    trail->UpdatePositions(std::move(positions));

    m_trail_node = make_unique<SceneNode>();
    m_trail_node->CreateComponent<ShadingComponent>(line_shader);
    auto& rendering_comp = m_trail_node->CreateComponent<RenderingComponent>(std::shared_ptr<VertexObject>(trail.release()));
    rendering_comp.SetDrawMode(DrawMode::Lines);
    AddChild(std::move(m_trail_node));
  }

  // Per-frame update
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

    // Update particle visual
    m_particle_node->GetTransform().SetPosition(m_state.positions[0]);

    // Update trail (add new position every few steps to avoid too many points)
    static int trail_counter = 0;
    if (++trail_counter % 10 == 0) {  // Add to trail every 10 steps
      m_trail_positions.push_back(m_state.positions[0]);

      // Limit trail length to prevent memory issues
      if (m_trail_positions.size() > 1000) {
        m_trail_positions.erase(m_trail_positions.begin());
      }

      // Update trail visualization
      auto positions = make_unique<PositionArray>(m_trail_positions);
      // Note: We would need to store a reference to the trail VertexObject to update it
      // For now, we'll just keep the trail static
    }
  }

  // Reset key
  void Reset() {
    m_state.positions[0] = glm::vec3(1.0f, 0.0f, 0.0f);
    m_state.velocities[0] = glm::vec3(0.0f, 0.0f, 0.0f);
    m_time = 0.0f;
    m_accumulator = 0.0;
    m_particle_node->GetTransform().SetPosition(m_state.positions[0]);

    // Reset trail
    m_trail_positions.clear();
    m_trail_positions.push_back(m_state.positions[0]);
  }

 private:
  float m_dt, m_time;
  double m_accumulator;

  ParticleState m_state;
  SimpleSystem m_system;

  std::unique_ptr<IntegratorBase<SimpleSystem, ParticleState>> m_integrator;
  SceneNode* m_particle_node;
  std::unique_ptr<SceneNode> m_trail_node;
  std::vector<glm::vec3> m_trail_positions;
};

}  // namespace GLOO
