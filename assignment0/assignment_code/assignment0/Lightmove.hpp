#pragma once
#include "gloo/SceneNode.hpp"

namespace GLOO {
class LightControlNode : public SceneNode {
 public:
  LightControlNode();
  void Update(double delta_time) override;

 private:
  float move_speed_;
};
}  // namespace GLOO
