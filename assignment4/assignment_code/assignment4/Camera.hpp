#ifndef CAMERA_H_
#define CAMERA_H_

#include "Ray.hpp"

namespace GLOO {
class Camera {
 public:
  virtual ~Camera() = default;
  virtual Ray GenerateRay(const glm::vec2& point) const = 0;
  virtual float GetTMin() const = 0;
};
}  // namespace GLOO

#endif
