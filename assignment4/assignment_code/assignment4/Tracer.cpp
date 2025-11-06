#include "Tracer.hpp"

#include <glm/gtx/string_cast.hpp>
#include <stdexcept>
#include <algorithm>
#include <limits>

#include "gloo/Transform.hpp"
#include "gloo/components/MaterialComponent.hpp"
#include "gloo/lights/AmbientLight.hpp"

#include "gloo/Image.hpp"
#include "Illuminator.hpp"

namespace GLOO {
void Tracer::Render(const Scene& scene, const std::string& output_file) {
  scene_ptr_ = &scene;

  auto& root = scene_ptr_->GetRootNode();
  tracing_components_ = root.GetComponentPtrsInChildren<TracingComponent>();
  light_components_ = root.GetComponentPtrsInChildren<LightComponent>();


  Image image(image_size_.x, image_size_.y);

  // Generate rays for each pixel and trace them
  for (size_t y = 0; y < image_size_.y; y++) {
    for (size_t x = 0; x < image_size_.x; x++) {
      // Convert pixel coordinates to normalized device coordinates [-1, 1]
      float u = 1.0f - (2.0f * x / image_size_.x);  // Flip x-axis
      float v = 1.0f - (2.0f * y / image_size_.y);  // Flip y-axis
      //flip y axis back
      v = -v;
      u = -u;

      // Generate ray from camera
      Ray ray = camera_.GenerateRay(glm::vec2(u, v));

      // Trace ray and get color (start with depth 0)
      HitRecord record;
      glm::vec3 color = TraceRay(ray, 0, record);

      // Set pixel color
      image.SetPixel(x, y, color);
    }
  }

  if (output_file.size())
    image.SavePNG(output_file);
}


glm::vec3 Tracer::TraceRay(const Ray& ray,
                           size_t bounces,
                           HitRecord& record) const {
  // Base case: if we've exceeded max recursion depth, return background color
  // Note: max_bounces_ = 0 means no recursion (just camera rays), so bounces=0 is still valid
  if (bounces > max_bounces_) {
    return GetBackgroundColor(ray.GetDirection());
  }

  // Initialize record to find closest hit
  record.time = std::numeric_limits<float>::max();
  TracingComponent* hit_component = nullptr;

  // Find closest intersection by testing all TracingComponents
  for (TracingComponent* tracing_comp : tracing_components_) {
    SceneNode* node = tracing_comp->GetNodePtr();
    const Transform& transform = node->GetTransform();

    // Transform ray from world space to object local space
    glm::mat4 world_to_local = glm::inverse(transform.GetLocalToWorldMatrix());

    // Transform ray origin (as point, w=1)
    glm::vec4 local_origin_homogeneous = world_to_local * glm::vec4(ray.GetOrigin(), 1.0f);
    glm::vec3 local_origin = glm::vec3(local_origin_homogeneous) / local_origin_homogeneous.w;

    // Transform ray direction (as vector, w=0) and normalize
    glm::vec4 local_dir_homogeneous = world_to_local * glm::vec4(ray.GetDirection(), 0.0f);
    glm::vec3 local_dir = glm::normalize(glm::vec3(local_dir_homogeneous));

    Ray local_ray(local_origin, local_dir);

    // Test intersection in local space
    HitRecord local_record;
    if (tracing_comp->GetHittable().Intersect(local_ray, camera_.GetTMin(), local_record)) {
      // Transform hit point back to world space
      glm::vec3 local_hit_point = local_ray.At(local_record.time);
      glm::vec4 world_hit_point_homogeneous = transform.GetLocalToWorldMatrix() * glm::vec4(local_hit_point, 1.0f);
      glm::vec3 world_hit_point = glm::vec3(world_hit_point_homogeneous) / world_hit_point_homogeneous.w;

      // Compute world-space t
      glm::vec3 ray_dir = ray.GetDirection();
      glm::vec3 to_hit = world_hit_point - ray.GetOrigin();
      float world_t = glm::dot(to_hit, ray_dir);

      // Ensure world_t is positive and closer than previous hits
      if (world_t > camera_.GetTMin() && world_t < record.time) {
        record.time = world_t;
        hit_component = tracing_comp;
        // Store local normal for now, will transform to world space later
        record.normal = local_record.normal;
      }
    }
  }

  // If no hit, return background color
  if (hit_component == nullptr) {
    return GetBackgroundColor(ray.GetDirection());
  }

  // Get material from the node that was hit
  SceneNode* hit_node = hit_component->GetNodePtr();
  MaterialComponent* material_comp = hit_node->GetComponentPtr<MaterialComponent>();
  const Material& material = material_comp ? material_comp->GetMaterial() : Material::GetDefault();

  // Transform normal to world space using inverse transpose
  const Transform& transform = hit_node->GetTransform();
  glm::mat4 local_to_world = transform.GetLocalToWorldMatrix();
  glm::mat4 inv_trans = glm::transpose(glm::inverse(local_to_world));
  glm::vec4 normal_homogeneous = inv_trans * glm::vec4(record.normal, 0.0f);
  glm::vec3 N = glm::normalize(glm::vec3(normal_homogeneous));

  // Compute hit point in world space
  glm::vec3 hit_point = ray.At(record.time);

  // Compute eye direction (direction from hit point back to the eye)
  glm::vec3 E = glm::normalize(-ray.GetDirection());

  // Initialize color with ambient term: I_ambient = L_ambient * k_ambient
  glm::vec3 ambient_light_color(0.0f);
  for (LightComponent* light_comp : light_components_) {
    if (light_comp->GetLightPtr()->GetType() == LightType::Ambient) {
      auto ambient_light_ptr = static_cast<AmbientLight*>(light_comp->GetLightPtr());
      ambient_light_color = ambient_light_ptr->GetAmbientColor();
      break;  // Typically only one ambient light
    }
  }
  glm::vec3 color = ambient_light_color * material.GetAmbientColor();

  // Accumulate diffuse and specular contributions from all non-ambient lights
  for (LightComponent* light_comp : light_components_) {
    // Skip ambient lights (already handled above)
    if (light_comp->GetLightPtr()->GetType() == LightType::Ambient) {
      continue;
    }

    // Get illumination information
    glm::vec3 L;
    glm::vec3 light_intensity;
    float dist_to_light;
    Illuminator::GetIllumination(*light_comp, hit_point, L, light_intensity, dist_to_light);

    // Shadow ray visibility testing: check if light is blocked
    bool light_visible = true;
    if (shadows_enabled_) {
      // Construct shadow ray: offset origin slightly to avoid self-intersection
      const float shadow_epsilon = 1e-4f;
      glm::vec3 shadow_origin = hit_point + shadow_epsilon * L;
      Ray shadow_ray(shadow_origin, L);

      // Test if any object intersects the shadow ray before reaching the light
      for (TracingComponent* tracing_comp : tracing_components_) {
        SceneNode* node = tracing_comp->GetNodePtr();
        const Transform& transform = node->GetTransform();

        // Transform shadow ray to object local space
        glm::mat4 world_to_local = glm::inverse(transform.GetLocalToWorldMatrix());

        glm::vec4 local_origin_homogeneous = world_to_local * glm::vec4(shadow_ray.GetOrigin(), 1.0f);
        glm::vec3 local_origin = glm::vec3(local_origin_homogeneous) / local_origin_homogeneous.w;

        glm::vec4 local_dir_homogeneous = world_to_local * glm::vec4(shadow_ray.GetDirection(), 0.0f);
        glm::vec3 local_dir = glm::normalize(glm::vec3(local_dir_homogeneous));

        Ray local_shadow_ray(local_origin, local_dir);

        // Test intersection in local space
        HitRecord shadow_record;
        if (tracing_comp->GetHittable().Intersect(local_shadow_ray, camera_.GetTMin(), shadow_record)) {
          // Transform hit point back to world space to compute world-space distance
          glm::vec3 local_hit_point = local_shadow_ray.At(shadow_record.time);
          glm::vec4 world_hit_point_homogeneous = transform.GetLocalToWorldMatrix() * glm::vec4(local_hit_point, 1.0f);
          glm::vec3 world_hit_point = glm::vec3(world_hit_point_homogeneous) / world_hit_point_homogeneous.w;

          // Compute distance along shadow ray
          glm::vec3 to_shadow_hit = world_hit_point - shadow_origin;
          float shadow_t = glm::dot(to_shadow_hit, L);

          // If intersection is before the light, light is blocked
          if (shadow_t > camera_.GetTMin() && shadow_t < dist_to_light) {
            light_visible = false;
            break;
          }
        }
      }
    }

    // Only add lighting contributions if light is visible (not in shadow)
    if (light_visible) {
      // Compute diffuse contribution: I_diffuse = max(N·L, 0) * lightIntensity * k_diffuse
      float ndotl = glm::dot(N, L);
      ndotl = std::max(ndotl, 0.0f);
      glm::vec3 I_diffuse = ndotl * light_intensity * material.GetDiffuseColor();

      // Compute specular contribution: I_spec = pow(max(R·L, 0), shininess) * lightIntensity * k_specular
      glm::vec3 R = glm::reflect(-E, N);  // Reflection of -E about N
      float rdotl = glm::dot(R, L);
      rdotl = std::max(rdotl, 0.0f);
      float shininess = material.GetShininess();
      glm::vec3 I_specular = std::pow(rdotl, shininess) * light_intensity * material.GetSpecularColor();

      // Add both diffuse and specular contributions
      color += I_diffuse + I_specular;
    }
  }

  // Compute direct lighting (ambient + diffuse + specular from lights)
  glm::vec3 I_direct = color;

  // Recursive specular reflection: if material has specular component, trace reflected ray
  glm::vec3 k_specular = material.GetSpecularColor();
  float specular_magnitude = glm::length(k_specular);

  if (specular_magnitude > 1e-6f && (bounces + 1) <= max_bounces_) {
    // Compute perfect reflection direction: R = reflect(-E, N)
    // Use the same eye direction E as in Phong specular for consistency
    // E points from hit point back to camera, so -E is the incident direction
    glm::vec3 R = glm::normalize(glm::reflect(-E, N));

    // Offset origin slightly along normal to avoid self-intersection
    const float offset_epsilon = 1e-4f;
    glm::vec3 new_origin = hit_point + offset_epsilon * R;

    // Spawn secondary reflected ray
    Ray reflected_ray(new_origin, R);

    // Recursively trace the reflected ray
    HitRecord reflected_record;
    glm::vec3 I_indirect = TraceRay(reflected_ray, bounces + 1, reflected_record);

    // Combine direct and indirect lighting: color = I_direct + k_specular * I_indirect
    color = I_direct + k_specular * I_indirect;
  } else {
    color = I_direct;
  }

  // Clamp color to [0, 1] per channel
  color = glm::clamp(color, 0.0f, 1.0f);

  return color;
}


glm::vec3 Tracer::GetBackgroundColor(const glm::vec3& direction) const {
  if (cube_map_ != nullptr) {
    return cube_map_->GetTexel(direction);
  } else
    return background_color_;
}
}  // namespace GLOO
