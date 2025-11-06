#include "Plane.hpp"

#include <cmath>
#include <glm/geometric.hpp>

namespace GLOO {
Plane::Plane(const glm::vec3& normal, float d) : normal_(glm::normalize(normal)), d_(d) {
}

bool Plane::Intersect(const Ray& ray, float t_min, HitRecord& record) const {
  // Ray equation: r(t) = ray.origin + t * ray.direction
  // Plane equation: n · P = d, where P is any point on the plane

  glm::vec3 ray_dir = ray.GetDirection();
  float ndotd = glm::dot(normal_, ray_dir);

  // Check if ray is parallel to plane (n · direction ≈ 0)
  const float epsilon = 1e-6f;
  if (std::abs(ndotd) < epsilon) {
    return false;  // Ray is parallel to plane, no intersection
  }

  // Solve for t: n · (origin + t * direction) = d
  // n · origin + t * (n · direction) = d
  // t = (d - n · origin) / (n · direction)
  float ndoto = glm::dot(normal_, ray.GetOrigin());
  float t = (d_ - ndoto) / ndotd;

  // Check if intersection is valid (t > t_min and closer than previous hit)
  if (t > t_min && t < record.time) {
    record.time = t;

    // Ensure normal faces opposite the incoming ray
    // If dot(n, direction) > 0, ray is coming from front, so flip normal
    if (ndotd > 0.0f) {
      record.normal = -normal_;
    } else {
      record.normal = normal_;
    }

    return true;
  }

  return false;
}
}  // namespace GLOO
