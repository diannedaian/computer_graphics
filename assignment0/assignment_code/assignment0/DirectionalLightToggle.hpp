#pragma once
#include "gloo/SceneNode.hpp"

namespace GLOO {
class DirectionalLightToggle : public SceneNode {
 public:
  DirectionalLightToggle();
  void Update(double delta_time) override;
  void SetDirectionalLightNode(SceneNode* directional_light_node);

 private:
  void ToggleDirectionalLight();
  SceneNode* directional_light_node_;
  bool is_enabled_;
};
}  // namespace GLOO
