#include "Triangle.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include <glm/common.hpp>
#include <glm/gtx/string_cast.hpp>

#include "Plane.hpp"

namespace GLOO {
Triangle::Triangle(const glm::vec3& p0,
                   const glm::vec3& p1,
                   const glm::vec3& p2,
                   const glm::vec3& n0,
                   const glm::vec3& n1,
                   const glm::vec3& n2) {
  positions_ = {p0, p1, p2};
  normals_ = {n0, n1, n2};
}

Triangle::Triangle(const std::vector<glm::vec3>& positions,
                   const std::vector<glm::vec3>& normals) {
  if (positions.size() != 3 || normals.size() != 3) {
    throw std::runtime_error("Triangle must have exactly 3 vertices and 3 normals");
  }
  positions_ = positions;
  normals_ = normals;
}

bool Triangle::Intersect(const Ray& ray, float t_min, HitRecord& record) const {
  // Ray-triangle intersection using barycentric coordinates
  // Solve A x = b, where A = [v1-v0, v2-v0, -ray.dir], x = [u, v, t], b = ray.origin - v0

  const glm::vec3& v0 = positions_[0];
  const glm::vec3& v1 = positions_[1];
  const glm::vec3& v2 = positions_[2];

  glm::vec3 edge1 = v1 - v0;
  glm::vec3 edge2 = v2 - v0;
  glm::vec3 ray_dir = ray.GetDirection();

  // Build matrix A = [edge1, edge2, -ray_dir]
  glm::mat3 A(edge1, edge2, -ray_dir);

  // Check if matrix is invertible (determinant not too small)
  const float epsilon = 1e-6f;
  float det = glm::determinant(A);
  if (std::abs(det) < epsilon) {
    return false;  // Ray is parallel to triangle plane
  }

  // Solve for [u, v, t] = A^(-1) * b, where b = ray.origin - v0
  glm::vec3 b = ray.GetOrigin() - v0;
  glm::mat3 A_inv = glm::inverse(A);
  glm::vec3 solution = A_inv * b;

  float u = solution.x;
  float v = solution.y;
  float t = solution.z;

  // Check if intersection is valid:
  // - u, v >= 0 (inside triangle)
  // - u + v <= 1 (inside triangle)
  // - t > t_min (forward along ray)
  // - t < record.time (closer than previous hit)
  if (u >= 0.0f && v >= 0.0f && (u + v) <= 1.0f && t > t_min && t < record.time) {
    record.time = t;

    // Interpolate normal: N = (1-u-v)*n0 + u*n1 + v*n2
    glm::vec3 n0 = normals_[0];
    glm::vec3 n1 = normals_[1];
    glm::vec3 n2 = normals_[2];
    float w = 1.0f - u - v;  // Barycentric weight for v0
    record.normal = glm::normalize(w * n0 + u * n1 + v * n2);

    return true;
  }

  return false;
}
}  // namespace GLOO
