#ifndef SPIRAL_CAMERA_H_
#define SPIRAL_CAMERA_H_

#include <cmath>
#include <glm/ext/quaternion_geometric.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "gloo/utils.hpp"

#include "Camera.hpp"
#include "Ray.hpp"
#include "CameraSpec.hpp"

namespace GLOO {
class SpiralCamera : public Camera {
 public:
  SpiralCamera(const glm::vec3& center,
               const glm::vec3& direction,
               const glm::vec3& up,
               float fov,
               float twistStrength);

  Ray GenerateRay(const glm::vec2& point) const override {
    // Extract x and y from point (normalized coordinates in [0, 1] from Tracer)
    float x = point.x;
    float y = point.y;

    // Compute normalized pixel offsets in [-1, 1]
    float nx = (2.0f * x - 1.0f);
    float ny = (2.0f * y - 1.0f);

    // Scale by FOV to account for field of view
    float d = 1.0f / tanf(fov_radian_ / 2.0f);

    // Compute a base direction (standard pinhole camera)
    // Note: w_ points forward (toward scene), so we use +d*w_ to get correct direction
    glm::vec3 dir = glm::normalize(nx * u_ + ny * v_ + d * w_);

    // Compute the rotation angle based on the horizontal coordinate
    // twistStrength is in radians, twist more near edges (nx = -1 or +1)
    float angle = twistStrength_ * nx;

    // Construct an orthonormal rotation matrix around w (camera forward axis)
    glm::mat3 rot = glm::mat3(glm::rotate(glm::mat4(1.0f), angle, w_));

    // Apply rotation to the direction
    dir = glm::normalize(rot * dir);

    // Return the ray
    return Ray(center_, dir);
  }

  float GetTMin() const override {
    return 0.0f;
  }

 private:
  float twistStrength_;  // radians of twist across the screen
  glm::vec3 center_;
  glm::vec3 u_, v_, w_;  // camera basis (right, up, forward)
  float fov_radian_;
};
}  // namespace GLOO

#endif
