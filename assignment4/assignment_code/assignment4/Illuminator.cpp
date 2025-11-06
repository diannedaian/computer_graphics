#include "Illuminator.hpp"

#include <limits>
#include <stdexcept>

#include <glm/geometric.hpp>

#include "gloo/lights/DirectionalLight.hpp"
#include "gloo/lights/PointLight.hpp"
#include "gloo/SceneNode.hpp"

namespace GLOO {
void Illuminator::GetIllumination(const LightComponent& light_component,
                                  const glm::vec3& hit_pos,
                                  glm::vec3& dir_to_light,
                                  glm::vec3& intensity,
                                  float& dist_to_light) {
  // Calculation will be done in world space.

  auto light_ptr = light_component.GetLightPtr();
  if (light_ptr->GetType() == LightType::Directional) {
    auto directional_light_ptr = static_cast<DirectionalLight*>(light_ptr);
    dir_to_light = -directional_light_ptr->GetDirection();
    intensity = directional_light_ptr->GetDiffuseColor();
    dist_to_light = std::numeric_limits<float>::max();
  } else if (light_ptr->GetType() == LightType::Point) {
    auto point_light_ptr = static_cast<PointLight*>(light_ptr);

    // Compute direction and distance from surface point to light position
    glm::vec3 light_pos = light_component.GetNodePtr()->GetTransform().GetWorldPosition();
    glm::vec3 to_light = light_pos - hit_pos;
    dist_to_light = glm::length(to_light);

    // Avoid divide-by-zero for very small distances
    const float epsilon = 1e-6f;
    if (dist_to_light < epsilon) {
      dist_to_light = epsilon;
    }
    dir_to_light = glm::normalize(to_light);

    // Compute intensity using inverse square law: I / (alpha * d^2)
    float alpha = point_light_ptr->GetAttenuation().x;
    glm::vec3 light_color = point_light_ptr->GetDiffuseColor();
    float d_squared = dist_to_light * dist_to_light;
    intensity = light_color / (alpha * d_squared);
  } else {
    throw std::runtime_error(
        "Unrecognized light type when computing "
        "illumination");
  }
}
}  // namespace GLOO
