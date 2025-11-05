#include "Plane.hpp"

namespace GLOO {
Plane::Plane(const glm::vec3& normal, float d) {
}

bool Plane::Intersect(const Ray& ray, float t_min, HitRecord& record) const {
  // TODO: Implement ray-plane intersection.
  return false;
}
}  // namespace GLOO
