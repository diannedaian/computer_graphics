#include "SpiralCamera.hpp"

namespace GLOO {
SpiralCamera::SpiralCamera(const glm::vec3& center,
                           const glm::vec3& direction,
                           const glm::vec3& up,
                           float fov,
                           float twistStrength)
    : twistStrength_(twistStrength),
      center_(center),
      fov_radian_(ToRadian(fov)) {
  // Normalize and set up the camera basis
  w_ = glm::normalize(direction);
  u_ = glm::normalize(glm::cross(w_, up));
  v_ = glm::normalize(glm::cross(u_, w_));
}
}  // namespace GLOO
