#include "Triangle.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include <glm/common.hpp>
#include <glm/gtx/string_cast.hpp>

#include "Plane.hpp"

namespace {
// Helper function to build an orthonormal basis from a normal vector
// Uses the standard "pick an arbitrary axis not parallel to N" method
void BuildOrthonormalBasisFromNormal(const glm::vec3& N, glm::vec3& T, glm::vec3& B) {
  // Normalize the input normal (should already be normalized, but be safe)
  glm::vec3 n = glm::normalize(N);

  // Choose an arbitrary axis that's not parallel to N
  // Use (1, 0, 0) if N is not too close to it, otherwise use (0, 1, 0)
  glm::vec3 axis = glm::abs(n.x) < 0.9f ? glm::vec3(1.0f, 0.0f, 0.0f) : glm::vec3(0.0f, 1.0f, 0.0f);

  // Compute tangent: T = normalize(axis - (axis · N) * N)
  T = glm::normalize(axis - glm::dot(axis, n) * n);

  // Compute bitangent: B = N × T
  B = glm::normalize(glm::cross(n, T));
}
}  // namespace

namespace GLOO {
Triangle::Triangle(const glm::vec3& p0,
                   const glm::vec3& p1,
                   const glm::vec3& p2,
                   const glm::vec3& n0,
                   const glm::vec3& n1,
                   const glm::vec3& n2) {
  positions_ = {p0, p1, p2};
  normals_ = {n0, n1, n2};
  uvs_ = {glm::vec2(0.0f, 0.0f), glm::vec2(1.0f, 0.0f), glm::vec2(0.5f, 1.0f)};  // Default UVs
  has_uvs_ = false;
}

Triangle::Triangle(const glm::vec3& p0,
                   const glm::vec3& p1,
                   const glm::vec3& p2,
                   const glm::vec3& n0,
                   const glm::vec3& n1,
                   const glm::vec3& n2,
                   const glm::vec2& uv0,
                   const glm::vec2& uv1,
                   const glm::vec2& uv2) {
  positions_ = {p0, p1, p2};
  normals_ = {n0, n1, n2};
  uvs_ = {uv0, uv1, uv2};
  has_uvs_ = true;
}

Triangle::Triangle(const std::vector<glm::vec3>& positions,
                   const std::vector<glm::vec3>& normals) {
  if (positions.size() != 3 || normals.size() != 3) {
    throw std::runtime_error("Triangle must have exactly 3 vertices and 3 normals");
  }
  positions_ = positions;
  normals_ = normals;
  uvs_ = {glm::vec2(0.0f, 0.0f), glm::vec2(1.0f, 0.0f), glm::vec2(0.5f, 1.0f)};  // Default UVs
  has_uvs_ = false;
}

Triangle::Triangle(const std::vector<glm::vec3>& positions,
                   const std::vector<glm::vec3>& normals,
                   const std::vector<glm::vec2>& uvs) {
  if (positions.size() != 3 || normals.size() != 3 || uvs.size() != 3) {
    throw std::runtime_error("Triangle must have exactly 3 vertices, 3 normals, and 3 UVs");
  }
  positions_ = positions;
  normals_ = normals;
  uvs_ = uvs;
  has_uvs_ = true;
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

    // Interpolate texture coordinates using barycentric coordinates
    // uv = w * uv0 + u * uv1 + v * uv2
    glm::vec2 uv0 = uvs_[0];
    glm::vec2 uv1 = uvs_[1];
    glm::vec2 uv2 = uvs_[2];
    record.uv = w * uv0 + u * uv1 + v * uv2;

    // Compute tangent frame (T, B, N) using UV derivatives
    glm::vec3 p0 = positions_[0];
    glm::vec3 p1 = positions_[1];
    glm::vec3 p2 = positions_[2];

    glm::vec3 dp1 = p1 - p0;
    glm::vec3 dp2 = p2 - p0;
    glm::vec2 duv1 = uv1 - uv0;
    glm::vec2 duv2 = uv2 - uv0;

    float det = duv1.x * duv2.y - duv1.y * duv2.x;
    const float epsilon = 1e-8f;

    if (std::abs(det) > epsilon) {
      float r = 1.0f / det;
      glm::vec3 T = (dp1 * duv2.y - dp2 * duv1.y) * r;
      glm::vec3 B = (dp2 * duv1.x - dp1 * duv2.x) * r;
      T = glm::normalize(T);
      B = glm::normalize(B);
      record.tangent = T;
      record.bitangent = B;
    } else {
      // Fallback: build any orthonormal frame from the normal
      BuildOrthonormalBasisFromNormal(record.normal, record.tangent, record.bitangent);
    }

    return true;
  }

  return false;
}
}  // namespace GLOO
