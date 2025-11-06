#include "Tracer.hpp"

#include <glm/gtx/string_cast.hpp>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <random>

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

  // Generate rays for each pixel with supersampling
  for (size_t y = 0; y < image_size_.y; y++) {
    for (size_t x = 0; x < image_size_.x; x++) {
      // Accumulate weighted color from all sub-samples
      glm::vec3 accum(0.0f);
      float wsum = 0.0f;

      // Loop over sub-samples within the pixel
      for (int sy = 0; sy < samples_per_dim_; ++sy) {
        for (int sx = 0; sx < samples_per_dim_; ++sx) {
          // Compute relative offset within pixel: dx, dy in [-0.5, 0.5]
          float dx = (static_cast<float>(sx) + 0.5f) / static_cast<float>(samples_per_dim_) - 0.5f;
          float dy = (static_cast<float>(sy) + 0.5f) / static_cast<float>(samples_per_dim_) - 0.5f;

          // Compute filter weight for this sample
          float w = ComputeFilterWeight(dx, dy, filter_type_);

          // Compute absolute pixel position for ray generation
          float pixel_x = static_cast<float>(x) + (static_cast<float>(sx) + 0.5f) / static_cast<float>(samples_per_dim_);
          float pixel_y = static_cast<float>(y) + (static_cast<float>(sy) + 0.5f) / static_cast<float>(samples_per_dim_);

          // Convert pixel coordinates to normalized device coordinates [-1, 1]
          // Normalize to [0, 1] first: pixel_x / width, pixel_y / height
          float u_norm = pixel_x / static_cast<float>(image_size_.x);
          float v_norm = pixel_y / static_cast<float>(image_size_.y);

          // Convert to [-1, 1] range and apply same flips as before
          float u = 1.0f - (2.0f * u_norm);  // Flip x-axis
          float v = 1.0f - (2.0f * v_norm);  // Flip y-axis
          // Flip back
          v = -v;
          u = -u;

          // Generate ray from camera
          Ray ray = camera_->GenerateRay(glm::vec2(u, v));

          // Trace ray and accumulate weighted color
          HitRecord record;
          glm::vec3 sample_color = TraceRay(ray, 0, record);
          accum += w * sample_color;
          wsum += w;
        }
      }

      // Normalize by total weight (weighted average)
      glm::vec3 color = (wsum > 0.0f) ? (accum / wsum) : glm::vec3(0.0f);

      // Set pixel color
      image.SetPixel(x, y, color);
    }
  }

  if (output_file.size())
    image.SavePNG(output_file);
}

void Tracer::SetSamplesPerDim(int samples_per_dim) {
  // Ensure samples_per_dim is at least 1
  if (samples_per_dim < 1) {
    samples_per_dim_ = 1;
  } else {
    samples_per_dim_ = samples_per_dim;
  }
}

void Tracer::SetFilterType(FilterType filter_type) {
  filter_type_ = filter_type;
}

namespace {
// Hash function to convert 3D position to a seed for RNG
uint32_t HashPosition(const glm::vec3& pos) {
  // Hash the position components
  uint32_t h = static_cast<uint32_t>(std::hash<float>{}(pos.x));
  h ^= static_cast<uint32_t>(std::hash<float>{}(pos.y)) + 0x9e3779b9 + (h << 6) + (h >> 2);
  h ^= static_cast<uint32_t>(std::hash<float>{}(pos.z)) + 0x9e3779b9 + (h << 6) + (h >> 2);
  return h;
}

// Build orthonormal basis from a single vector (for R)
void BuildOrthonormalBasisFromVector(const glm::vec3& R, glm::vec3& T, glm::vec3& B) {
  glm::vec3 r = glm::normalize(R);

  // Choose an arbitrary axis that's not parallel to R
  glm::vec3 axis = glm::abs(r.x) < 0.9f ? glm::vec3(1.0f, 0.0f, 0.0f) : glm::vec3(0.0f, 1.0f, 0.0f);

  // Compute tangent: T = normalize(axis - (axis · R) * R)
  T = glm::normalize(axis - glm::dot(axis, r) * r);

  // Compute bitangent: B = R × T
  B = glm::normalize(glm::cross(r, T));
}

// Generate a random direction in spherical coordinates around R using seeded RNG
glm::vec3 SampleGlossyDirection(const glm::vec3& R, float roughness, const glm::vec3& hit_point, int sample_index) {
  // Create seeded RNG based on hit point and sample index for deterministic per-pixel sampling
  uint32_t seed = HashPosition(hit_point) + static_cast<uint32_t>(sample_index);
  std::mt19937 gen(seed);
  std::uniform_real_distribution<float> dis(0.0f, 1.0f);

  // Generate random values for spherical coordinates
  float rand1 = dis(gen);
  float rand2 = dis(gen);

  // Spherical coordinate jittering: theta based on roughness
  // Using formula: theta = acos(pow(rand1, 1.0f/(roughness+1e-3f)))
  float theta = std::acos(std::pow(rand1, 1.0f / (roughness + 1e-3f)));
  float phi = 2.0f * M_PI * rand2;

  // Build orthonormal basis around R (T, B, R)
  glm::vec3 T, B;
  BuildOrthonormalBasisFromVector(R, T, B);

  // Convert spherical coordinates to vector in local frame
  // In local frame: x = sin(theta)*cos(phi), y = sin(theta)*sin(phi), z = cos(theta)
  // Where z-axis is R, x-axis is T, y-axis is B
  glm::vec3 local_dir(
    std::sin(theta) * std::cos(phi),  // T component
    std::sin(theta) * std::sin(phi),  // B component
    std::cos(theta)                   // R component
  );

  // Transform from local frame to world space using TBN basis
  glm::vec3 world_dir = local_dir.x * T + local_dir.y * B + local_dir.z * R;

  return glm::normalize(world_dir);
}
}  // namespace

float Tracer::ComputeFilterWeight(float dx, float dy, FilterType filter_type) const {
  switch (filter_type) {
    case FilterType::Box:
      // Box filter: uniform weight for all samples
      return 1.0f;

    case FilterType::Tent: {
      // Tent filter: higher weight near pixel center
      // Weight decreases linearly with distance from center
      float wx = std::max(0.0f, 1.0f - std::abs(dx) * 2.0f);
      float wy = std::max(0.0f, 1.0f - std::abs(dy) * 2.0f);
      return wx * wy;
    }

    default:
      return 1.0f;  // Default to box filter
  }
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
        if (tracing_comp->GetHittable().Intersect(local_ray, camera_->GetTMin(), local_record)) {
      // Transform hit point back to world space
      glm::vec3 local_hit_point = local_ray.At(local_record.time);
      glm::vec4 world_hit_point_homogeneous = transform.GetLocalToWorldMatrix() * glm::vec4(local_hit_point, 1.0f);
      glm::vec3 world_hit_point = glm::vec3(world_hit_point_homogeneous) / world_hit_point_homogeneous.w;

      // Compute world-space t
      glm::vec3 ray_dir = ray.GetDirection();
      glm::vec3 to_hit = world_hit_point - ray.GetOrigin();
      float world_t = glm::dot(to_hit, ray_dir);

      // Ensure world_t is positive and closer than previous hits
      if (world_t > camera_->GetTMin() && world_t < record.time) {
        record.time = world_t;
        hit_component = tracing_comp;
        // Store local hit record fields (normal, tangent, bitangent, uv) for later transformation
        record.normal = local_record.normal;
        record.tangent = local_record.tangent;
        record.bitangent = local_record.bitangent;
        record.uv = local_record.uv;
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

  // Transform normal to world space
  glm::vec4 normal_homogeneous = inv_trans * glm::vec4(record.normal, 0.0f);
  glm::vec3 N = glm::normalize(glm::vec3(normal_homogeneous));

  // Apply normal mapping if material has a normal map
  if (material.HasNormalMap()) {
    // Only transform tangent and bitangent when we need them for normal mapping
    // Transform tangent to world space
    glm::vec4 tangent_homogeneous = inv_trans * glm::vec4(record.tangent, 0.0f);
    glm::vec3 T = glm::normalize(glm::vec3(tangent_homogeneous));

    // Transform bitangent to world space
    glm::vec4 bitangent_homogeneous = inv_trans * glm::vec4(record.bitangent, 0.0f);
    glm::vec3 B = glm::normalize(glm::vec3(bitangent_homogeneous));
    // Sample tangent-space normal from normal map
    glm::vec3 n_ts = material.SampleNormalMap(record.uv.x, record.uv.y);

    // Build TBN matrix: [T, B, N] to transform from tangent space to world space
    // Ensure all vectors are normalized
    T = glm::normalize(T);
    B = glm::normalize(B);
    N = glm::normalize(N);

    // Check for degenerate cases (NaN or zero vectors)
    bool has_nan = std::isnan(T.x) || std::isnan(T.y) || std::isnan(T.z) ||
                   std::isnan(B.x) || std::isnan(B.y) || std::isnan(B.z) ||
                   std::isnan(N.x) || std::isnan(N.y) || std::isnan(N.z);
    bool has_zero = glm::length(T) < 1e-6f || glm::length(B) < 1e-6f || glm::length(N) < 1e-6f;

    if (!has_nan && !has_zero && glm::length(n_ts) > 1e-6f) {
      glm::mat3 TBN = glm::mat3(T, B, N);
      glm::vec3 N_perturbed = TBN * n_ts;
      float len = glm::length(N_perturbed);
      if (len > 1e-6f && !std::isnan(len)) {
        N = glm::normalize(N_perturbed);
      }
      // If perturbation failed, fall back to original N
    }
    // If normal mapping failed, use original geometric normal N
  }

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
        if (tracing_comp->GetHittable().Intersect(local_shadow_ray, camera_->GetTMin(), shadow_record)) {
          // Transform hit point back to world space to compute world-space distance
          glm::vec3 local_hit_point = local_shadow_ray.At(shadow_record.time);
          glm::vec4 world_hit_point_homogeneous = transform.GetLocalToWorldMatrix() * glm::vec4(local_hit_point, 1.0f);
          glm::vec3 world_hit_point = glm::vec3(world_hit_point_homogeneous) / world_hit_point_homogeneous.w;

          // Compute distance along shadow ray
          glm::vec3 to_shadow_hit = world_hit_point - shadow_origin;
          float shadow_t = glm::dot(to_shadow_hit, L);

          // If intersection is before the light, light is blocked
          if (shadow_t > camera_->GetTMin() && shadow_t < dist_to_light) {
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
    glm::vec3 new_origin = hit_point + offset_epsilon * N;

    glm::vec3 I_indirect;

    // Check if material is glossy (roughness > 0)
    if (material.IsGlossy()) {
      // Glossy reflection: sample multiple directions around R using Monte Carlo
      // Use 8-16 samples for reasonable performance
      const int glossy_samples = 8;
      float roughness = material.GetRoughness();

      glm::vec3 acc(0.0f);
      for (int i = 0; i < glossy_samples; i++) {
        // Generate random direction using spherical coordinates with seeded RNG
        // This ensures per-pixel deterministic sampling (same hit point = same samples)
        glm::vec3 random_dir = SampleGlossyDirection(R, roughness, hit_point, i);

        // Spawn reflected ray with random direction
        Ray reflected_ray(new_origin, random_dir);

        // Recursively trace the reflected ray
        HitRecord reflected_record;
        acc += TraceRay(reflected_ray, bounces + 1, reflected_record);
      }

      // Average the samples
      I_indirect = acc / static_cast<float>(glossy_samples);
    } else {
      // Perfect mirror reflection (roughness == 0)
      Ray reflected_ray(new_origin, R);

      // Recursively trace the reflected ray
      HitRecord reflected_record;
      I_indirect = TraceRay(reflected_ray, bounces + 1, reflected_record);
    }

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
