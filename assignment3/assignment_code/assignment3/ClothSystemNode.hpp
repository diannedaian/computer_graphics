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

class ClothSystemNode : public SceneNode {
 public:
  ClothSystemNode(IntegratorType integrator_type, float dt, int cloth_size, float rest_length)
      : m_dt(dt), m_time(0.0f), m_cloth_size(cloth_size), m_rest_length(rest_length), m_accumulator(0.0) {
    std::cerr << "Initializing ClothSystemNode with " << cloth_size << "x" << cloth_size << " cloth..." << std::endl;

    // ---- 1. Build the cloth system ----
    const float k_structural = 50.0f;   // Spring constant for structural springs
    const float k_shear = 25.0f;        // Spring constant for shear springs
    const float k_flex = 10.0f;         // Spring constant for flex springs

    // Add particles in a grid pattern
    for (int i = 0; i < m_cloth_size; i++) {
      for (int j = 0; j < m_cloth_size; j++) {
        glm::vec3 pos(j * m_rest_length, -i * m_rest_length, 0.0f);
        m_system.AddParticle(pos, 1.0f, glm::vec3(0.0f));
      }
    }

    // Add structural springs (horizontal and vertical connections)
    // Only add each spring once by iterating through each edge
    for (int i = 0; i < m_cloth_size; i++) {
      for (int j = 0; j < m_cloth_size; j++) {
        int current = IndexOf(i, j);

        // Horizontal structural springs (right neighbor) - only add from left to right
        if (j < m_cloth_size - 1) {
          int right = IndexOf(i, j + 1);
          m_system.AddSpring(current, right, k_structural, m_rest_length);
        }

        // Vertical structural springs (bottom neighbor) - only add from top to bottom
        if (i < m_cloth_size - 1) {
          int bottom = IndexOf(i + 1, j);
          m_system.AddSpring(current, bottom, k_structural, m_rest_length);
        }
      }
    }

    // Add shear springs (diagonal connections) - only add from top-left to bottom-right
    for (int i = 0; i < m_cloth_size - 1; i++) {
      for (int j = 0; j < m_cloth_size - 1; j++) {
        int current = IndexOf(i, j);
        int diag = IndexOf(i + 1, j + 1);
        m_system.AddSpring(current, diag, k_shear, m_rest_length * 1.414f); // sqrt(2) * rest_length
      }
    }

    // Add shear springs (diagonal connections) - only add from top-right to bottom-left
    for (int i = 0; i < m_cloth_size - 1; i++) {
      for (int j = 1; j < m_cloth_size; j++) {
        int current = IndexOf(i, j);
        int diag = IndexOf(i + 1, j - 1);
        m_system.AddSpring(current, diag, k_shear, m_rest_length * 1.414f);
      }
    }

    // Add flex springs (skip-one connections for bending resistance)
    for (int i = 0; i < m_cloth_size; i++) {
      for (int j = 0; j < m_cloth_size; j++) {
        int current = IndexOf(i, j);

        // Horizontal flex springs (skip one to the right)
        if (j < m_cloth_size - 2) {
          int right2 = IndexOf(i, j + 2);
          m_system.AddSpring(current, right2, k_flex, m_rest_length * 2.0f);
        }

        // Vertical flex springs (skip one down)
        if (i < m_cloth_size - 2) {
          int bottom2 = IndexOf(i + 2, j);
          m_system.AddSpring(current, bottom2, k_flex, m_rest_length * 2.0f);
        }
      }
    }

    // Fix the top row of particles
    for (int j = 0; j < m_cloth_size; j++) {
      m_system.FixParticle(IndexOf(0, j));
    }

    // Keep a copy of the starting state
    m_state = m_system.GetInitialState();

    // ---- 2. Create the integrator ----
    m_integrator = IntegratorFactory::CreateIntegrator<PendulumSystem, ParticleState>(integrator_type);

    // ---- 3. Build visual spheres for particles ----
    auto shader = std::make_shared<SimpleShader>();
    auto sphere = PrimitiveFactory::CreateSphere(0.02f, 15, 15);
    auto sphere_shared = std::shared_ptr<VertexObject>(sphere.release());

    for (int i = 0; i < m_cloth_size * m_cloth_size; i++) {
      auto child = make_unique<SceneNode>();
      child->CreateComponent<ShadingComponent>(shader);
      child->CreateComponent<RenderingComponent>(sphere_shared);
      child->GetTransform().SetPosition(m_state.positions[i]);
      SceneNode* ptr = child.get();
      AddChild(std::move(child));
      m_particle_nodes.push_back(ptr);
    }

    // ---- 4. Build visual wireframe for cloth ----
    CreateClothWireframe();

    // ---- 5. Create wind visualization ----
    CreateWindVisualization();

    std::cerr << "ClothSystemNode initialized with "
              << m_cloth_size * m_cloth_size << " particles and "
              << m_system.springs.size() << " springs." << std::endl;
  }

  // Helper method to convert 2D grid coordinates to 1D particle index
  int IndexOf(int i, int j) const {
    return i * m_cloth_size + j;
  }

  // Create wireframe visualization for the cloth
  void CreateClothWireframe() {
    auto line_shader = std::make_shared<SimpleShader>();

    // Build line segments for structural springs only (for cleaner wireframe)
    m_wireframe_positions.clear();

    // Add horizontal structural lines
    for (int i = 0; i < m_cloth_size; i++) {
      for (int j = 0; j < m_cloth_size - 1; j++) {
        m_wireframe_positions.push_back(m_state.positions[IndexOf(i, j)]);
        m_wireframe_positions.push_back(m_state.positions[IndexOf(i, j + 1)]);
      }
    }

    // Add vertical structural lines
    for (int i = 0; i < m_cloth_size - 1; i++) {
      for (int j = 0; j < m_cloth_size; j++) {
        m_wireframe_positions.push_back(m_state.positions[IndexOf(i, j)]);
        m_wireframe_positions.push_back(m_state.positions[IndexOf(i + 1, j)]);
      }
    }

    // Create VertexObject with the line data using the correct API
    auto wireframe = make_unique<VertexObject>();
    auto positions = make_unique<PositionArray>(m_wireframe_positions);
    wireframe->UpdatePositions(std::move(positions));

    // Store reference to the wireframe for updates
    m_wireframe_vertex_object = wireframe.get();

    // Create the wireframe node
    m_wireframe_node = make_unique<SceneNode>();
    m_wireframe_node->CreateComponent<ShadingComponent>(line_shader);
    auto& rendering_comp = m_wireframe_node->CreateComponent<RenderingComponent>(std::shared_ptr<VertexObject>(wireframe.release()));
    rendering_comp.SetDrawMode(DrawMode::Lines);
    AddChild(std::move(m_wireframe_node));
  }

  // Create wind visualization (arrow showing wind direction)
  void CreateWindVisualization() {
    auto line_shader = std::make_shared<SimpleShader>();

    // Create wind arrow (simple line with arrowhead)
    std::vector<glm::vec3> wind_arrow_positions;

    // Arrow shaft
    wind_arrow_positions.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
    wind_arrow_positions.push_back(glm::vec3(0.5f, 0.0f, 0.0f));

    // Arrow head
    wind_arrow_positions.push_back(glm::vec3(0.5f, 0.0f, 0.0f));
    wind_arrow_positions.push_back(glm::vec3(0.4f, 0.1f, 0.0f));

    wind_arrow_positions.push_back(glm::vec3(0.5f, 0.0f, 0.0f));
    wind_arrow_positions.push_back(glm::vec3(0.4f, -0.1f, 0.0f));

    // Create wind arrow object
    auto wind_arrow = make_unique<VertexObject>();
    auto positions = make_unique<PositionArray>(wind_arrow_positions);
    wind_arrow->UpdatePositions(std::move(positions));

    // Store reference for updates
    m_wind_arrow_object = wind_arrow.get();

    // Create wind visualization node
    m_wind_visualization_node = make_unique<SceneNode>();
    m_wind_visualization_node->CreateComponent<ShadingComponent>(line_shader);
    auto& rendering_comp = m_wind_visualization_node->CreateComponent<RenderingComponent>(std::shared_ptr<VertexObject>(wind_arrow.release()));
    rendering_comp.SetDrawMode(DrawMode::Lines);

    // Position wind arrow above the cloth
    m_wind_visualization_node->GetTransform().SetPosition(glm::vec3(0.0f, 1.0f, 0.0f));
    AddChild(std::move(m_wind_visualization_node));
  }

  // Update wireframe with current particle positions
  void UpdateWireframe() {
    if (!m_wireframe_vertex_object) return;

    // Update wireframe positions with current particle positions
    unsigned int index = 0;

    // Update horizontal structural lines
    for (int i = 0; i < m_cloth_size; i++) {
      for (int j = 0; j < m_cloth_size - 1; j++) {
        m_wireframe_positions[index++] = m_state.positions[IndexOf(i, j)];
        m_wireframe_positions[index++] = m_state.positions[IndexOf(i, j + 1)];
      }
    }

    // Update vertical structural lines
    for (int i = 0; i < m_cloth_size - 1; i++) {
      for (int j = 0; j < m_cloth_size; j++) {
        m_wireframe_positions[index++] = m_state.positions[IndexOf(i, j)];
        m_wireframe_positions[index++] = m_state.positions[IndexOf(i + 1, j)];
      }
    }

    // Update the vertex buffer using the correct API
    auto positions = make_unique<PositionArray>(m_wireframe_positions);
    m_wireframe_vertex_object->UpdatePositions(std::move(positions));
  }

  // Update wind visualization arrow
  void UpdateWindVisualization() {
    if (!m_wind_arrow_object || !m_wind_visualization_node) return;

    // Only show wind arrow when wind is enabled
    m_wind_visualization_node->GetTransform().SetPosition(
      m_system.IsWindEnabled() ? glm::vec3(0.0f, 1.0f, 0.0f) : glm::vec3(0.0f, -10.0f, 0.0f)
    );

    if (!m_system.IsWindEnabled()) return;

    // Update arrow direction based on wind direction
    glm::vec3 wind_dir = m_system.wind_direction;
    float wind_strength = m_system.wind_strength;

    // Create arrow pointing in wind direction
    std::vector<glm::vec3> wind_arrow_positions;

    // Arrow shaft
    wind_arrow_positions.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
    wind_arrow_positions.push_back(wind_dir * 0.5f * wind_strength);

    // Arrow head
    glm::vec3 arrow_tip = wind_dir * 0.5f * wind_strength;
    glm::vec3 perp = glm::normalize(glm::cross(wind_dir, glm::vec3(0.0f, 1.0f, 0.0f))) * 0.1f;

    wind_arrow_positions.push_back(arrow_tip);
    wind_arrow_positions.push_back(arrow_tip - wind_dir * 0.1f + perp);

    wind_arrow_positions.push_back(arrow_tip);
    wind_arrow_positions.push_back(arrow_tip - wind_dir * 0.1f - perp);

    // Update the wind arrow
    auto positions = make_unique<PositionArray>(wind_arrow_positions);
    m_wind_arrow_object->UpdatePositions(std::move(positions));
  }

  // Per-frame update
  void Update(double delta_time) override {
    if (InputManager::GetInstance().IsKeyPressed('R'))
      Reset();

    // Toggle wind with 'W' key
    if (InputManager::GetInstance().IsKeyPressed('W')) {
      m_system.ToggleWind();
      std::cerr << "Wind " << (m_system.IsWindEnabled() ? "enabled" : "disabled") << std::endl;
    }

    // Update wind
    m_system.UpdateWind(m_dt);

    // Accumulate time so simulation is frame-rate independent
    m_accumulator += delta_time;
    while (m_accumulator >= m_dt) {
      m_state = m_integrator->Integrate(m_system, m_state, m_time, m_dt);
      m_time += m_dt;
      m_accumulator -= m_dt;
    }

    // Update particle visuals
    for (int i = 0; i < m_cloth_size * m_cloth_size; i++) {
      m_particle_nodes[i]->GetTransform().SetPosition(m_state.positions[i]);
    }

    // Update wireframe (we'll recreate it each frame for simplicity)
    UpdateWireframe();

    // Update wind visualization
    UpdateWindVisualization();
  }

  // Reset key
  void Reset() {
    m_state = m_system.GetInitialState();
    m_time = 0.0f;
    m_accumulator = 0.0;
    for (int i = 0; i < m_cloth_size * m_cloth_size; i++) {
      m_particle_nodes[i]->GetTransform().SetPosition(m_state.positions[i]);
    }
  }

 private:
  float m_dt, m_time;
  double m_accumulator;
  int m_cloth_size;
  float m_rest_length;

  ParticleState m_state;
  PendulumSystem m_system;

  std::unique_ptr<IntegratorBase<PendulumSystem, ParticleState>> m_integrator;
  std::vector<SceneNode*> m_particle_nodes;
  std::unique_ptr<SceneNode> m_wireframe_node;
  VertexObject* m_wireframe_vertex_object;
  std::vector<glm::vec3> m_wireframe_positions;

  // Wind visualization
  std::unique_ptr<SceneNode> m_wind_visualization_node;
  VertexObject* m_wind_arrow_object;
};

}  // namespace GLOO
