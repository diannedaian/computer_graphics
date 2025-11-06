#ifndef GLOO_MATERIAL_H_
#define GLOO_MATERIAL_H_

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <memory>
#include <algorithm>
#include "gloo/Image.hpp"

namespace GLOO {
class Material {
 public:
  Material()
      : ambient_color_(0.0f),
        diffuse_color_(0.0f),
        specular_color_(0.0f),
        shininess_(0.0f),
        roughness_(0.0f) {
  }
  Material(const glm::vec3& ambient_color,
           const glm::vec3& diffuse_color,
           const glm::vec3& specular_color,
           float shininess)
      : ambient_color_(ambient_color),
        diffuse_color_(diffuse_color),
        specular_color_(specular_color),
        shininess_(shininess),
        roughness_(0.0f) {
  }

  static const Material& GetDefault() {
    static Material default_material(glm::vec3(0.5f, 0.1f, 0.2f),
                                     glm::vec3(0.5f, 0.1f, 0.2f),
                                     glm::vec3(0.4f, 0.4f, 0.4f), 20.0f);
    return default_material;
  }

  glm::vec3 GetAmbientColor() const {
    return ambient_color_;
  }

  void SetAmbientColor(const glm::vec3& color) {
    ambient_color_ = color;
  }

  glm::vec3 GetDiffuseColor() const {
    return diffuse_color_;
  }

  void SetDiffuseColor(const glm::vec3& color) {
    diffuse_color_ = color;
  }

  glm::vec3 GetSpecularColor() const {
    return specular_color_;
  }

  void SetSpecularColor(const glm::vec3& color) {
    specular_color_ = color;
  }

  float GetShininess() const {
    return shininess_;
  }

  void SetShininess(float shininess) {
    shininess_ = shininess;
  }

  // Normal map support
  bool HasNormalMap() const {
    return use_normal_map_ && normal_map_ != nullptr;
  }

  glm::vec3 SampleNormalMap(float u, float v) const {
    if (!HasNormalMap()) {
      return glm::vec3(0.0f, 0.0f, 1.0f);  // Default normal (pointing up)
    }

    // Clamp or wrap u, v into [0, 1]
    u = glm::clamp(u, 0.0f, 1.0f);
    v = glm::clamp(v, 0.0f, 1.0f);

    // Convert to pixel coordinates
    size_t width = normal_map_->GetWidth();
    size_t height = normal_map_->GetHeight();
    size_t ix = static_cast<size_t>(u * static_cast<float>(width));
    size_t iy = static_cast<size_t>(v * static_cast<float>(height));

    // Clamp to valid pixel range
    ix = std::min(ix, width - 1);
    iy = std::min(iy, height - 1);

    // Read the texel from normal map
    const glm::vec3& c = normal_map_->GetPixel(ix, iy);

    // Map RGB from [0, 1] to [-1, 1]
    glm::vec3 n_ts = glm::vec3(2.0f * c.r - 1.0f,
                               2.0f * c.g - 1.0f,
                               2.0f * c.b - 1.0f);

    return glm::normalize(n_ts);
  }

  void SetNormalMap(std::shared_ptr<Image> normal_map) {
    normal_map_ = normal_map;
  }

  void SetUseNormalMap(bool use_normal_map) {
    use_normal_map_ = use_normal_map;
  }

  // Roughness for glossy reflection (0 = perfect mirror, >0 = glossy)
  float GetRoughness() const {
    return roughness_;
  }

  void SetRoughness(float roughness) {
    roughness_ = roughness;
  }

  // Check if material is glossy (reflective with roughness > 0)
  bool IsGlossy() const {
    float specular_magnitude = glm::length(specular_color_);
    return specular_magnitude > 1e-6f && roughness_ > 0.0f;
  }

 private:
  glm::vec3 ambient_color_;
  glm::vec3 diffuse_color_;
  glm::vec3 specular_color_;
  float shininess_;

  std::shared_ptr<Image> normal_map_;
  bool use_normal_map_ = false;
  float roughness_ = 0.0f;  // 0 = perfect mirror, >0 = glossy
};
}  // namespace GLOO

#endif
