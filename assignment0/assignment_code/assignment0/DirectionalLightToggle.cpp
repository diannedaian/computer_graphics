#include "DirectionalLightToggle.hpp"
#include "gloo/InputManager.hpp"
#include "gloo/components/LightComponent.hpp"
#include "gloo/lights/DirectionalLight.hpp"

namespace GLOO {
DirectionalLightToggle::DirectionalLightToggle()
    : directional_light_node_(nullptr), is_enabled_(true) {
}

void DirectionalLightToggle::SetDirectionalLightNode(SceneNode* directional_light_node) {
  directional_light_node_ = directional_light_node;
}

void DirectionalLightToggle::Update(double delta_time) {
  // Handle 'P' key input for directional light toggling
  static bool prev_released = true;
  if (InputManager::GetInstance().IsKeyPressed('P')) {
    if (prev_released) {
      ToggleDirectionalLight();
    }
    prev_released = false;
  } else if (InputManager::GetInstance().IsKeyReleased('P')) {
    prev_released = true;
  }
}

void DirectionalLightToggle::ToggleDirectionalLight() {
  if (directional_light_node_ != nullptr) {
    is_enabled_ = !is_enabled_;
    // Toggle the light component by removing/adding it
    if (is_enabled_) {
      // Re-add the light component
      auto directional_light = std::make_shared<DirectionalLight>();
      directional_light->SetDirection(glm::vec3(0.0f, -1.0f, -0.5f));
      directional_light->SetDiffuseColor(glm::vec3(0.8f, 0.8f, 0.8f));
      directional_light->SetSpecularColor(glm::vec3(1.0f, 1.0f, 1.0f));
      directional_light_node_->CreateComponent<LightComponent>(directional_light);
    } else {
      // Remove the light component
      directional_light_node_->RemoveComponent<LightComponent>();
    }
  }
}
}  // namespace GLOO
