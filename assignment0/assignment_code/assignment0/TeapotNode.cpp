#include "TeapotNode.hpp"
#include "gloo/shaders/PhongShader.hpp"
#include "gloo/components/RenderingComponent.hpp"
#include "gloo/components/ShadingComponent.hpp"
#include "gloo/components/MaterialComponent.hpp"
#include "gloo/Material.hpp"
#include "gloo/MeshLoader.hpp"
#include "gloo/InputManager.hpp"

namespace GLOO {
TeapotNode::TeapotNode() {
  // Create and set up the Phong shader
  std::shared_ptr<PhongShader> shader = std::make_shared<PhongShader>();

  // Load the teapot mesh
  std::shared_ptr<VertexObject> mesh =
      MeshLoader::Import("assignment0/teapot.obj").vertex_obj;
  if (mesh == nullptr) {
    return;
  }

  // Create the shading, rendering, and material components
  CreateComponent<ShadingComponent>(shader);
  CreateComponent<RenderingComponent>(mesh);
  CreateComponent<MaterialComponent>(
      std::make_shared<Material>(Material::GetDefault()));

  // Set the teapot's position and rotation
  GetTransform().SetPosition(glm::vec3(0.f, 0.f, 0.f));
  GetTransform().SetRotation(glm::vec3(1.0f, 0.0f, 0.0f), 0.3f);
}

void TeapotNode::Update(double delta_time) {
  // Handle 'C' key input for color toggling
  static bool prev_released = true;
  if (InputManager::GetInstance().IsKeyPressed('C')) {
    if (prev_released) {
      ToggleColor();
    }
    prev_released = false;
  } else if (InputManager::GetInstance().IsKeyReleased('C')) {
    prev_released = true;
  }
}

void TeapotNode::ToggleColor() {
  // Get the material component and toggle between red and blue
  auto material_comp = GetComponentPtr<MaterialComponent>();
  if (material_comp != nullptr) {
    auto& material = material_comp->GetMaterial();

    // Toggle between red and blue
    static bool is_red = true; // Start with red
    if (is_red) {
      // Switch to blue
      material.SetAmbientColor(glm::vec3(0.0f, 0.0f, 0.2f));
      material.SetDiffuseColor(glm::vec3(0.0f, 0.0f, 0.8f));
      material.SetSpecularColor(glm::vec3(0.0f, 0.0f, 1.0f));
    } else {
      // Switch to red
      material.SetAmbientColor(glm::vec3(0.2f, 0.0f, 0.0f));
      material.SetDiffuseColor(glm::vec3(0.8f, 0.0f, 0.0f));
      material.SetSpecularColor(glm::vec3(1.0f, 0.0f, 0.0f));
    }
    is_red = !is_red;
  }
}
}
