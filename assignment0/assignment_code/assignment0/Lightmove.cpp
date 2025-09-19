#include "Lightmove.hpp"
#include "gloo/InputManager.hpp"
#include "gloo/external.hpp"

namespace GLOO {
LightControlNode::LightControlNode() : move_speed_(5.0f) {
}

void LightControlNode::Update(double delta_time) {
  float delta_dist = move_speed_ * static_cast<float>(delta_time);
  glm::vec3 current_pos = GetTransform().GetPosition();
  glm::vec3 new_pos = current_pos;

  // Handle arrow key input
  if (InputManager::GetInstance().IsKeyPressed(GLFW_KEY_UP)) {
    new_pos.y += delta_dist;
  }
  if (InputManager::GetInstance().IsKeyPressed(GLFW_KEY_DOWN)) {
    new_pos.y -= delta_dist;
  }
  if (InputManager::GetInstance().IsKeyPressed(GLFW_KEY_LEFT)) {
    new_pos.x -= delta_dist;
  }
  if (InputManager::GetInstance().IsKeyPressed(GLFW_KEY_RIGHT)) {
    new_pos.x += delta_dist;
  }

  GetTransform().SetPosition(new_pos);
}
}  // namespace GLOO
