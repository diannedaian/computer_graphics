#ifndef HIT_RECORD_H_
#define HIT_RECORD_H_

#include <limits>
#include <ostream>

#include <glm/gtx/string_cast.hpp>

#include "gloo/Material.hpp"

namespace GLOO {
struct HitRecord {
  HitRecord() {
    time = std::numeric_limits<float>::max();
    uv = glm::vec2(0.0f, 0.0f);  // Default UV coordinates
    tangent = glm::vec3(1.0f, 0.0f, 0.0f);   // Default tangent
    bitangent = glm::vec3(0.0f, 1.0f, 0.0f); // Default bitangent
  }

  float time;
  glm::vec3 normal;
  glm::vec2 uv;  // Texture coordinates for normal mapping
  glm::vec3 tangent;   // Tangent vector (T in TBN frame)
  glm::vec3 bitangent; // Bitangent vector (B in TBN frame)
};

inline std::ostream& operator<<(std::ostream& os, const HitRecord& rec) {
  os << "HitRecord <" << rec.time << ", " << glm::to_string(rec.normal) << ">";
  return os;
}

}  // namespace GLOO

#endif
