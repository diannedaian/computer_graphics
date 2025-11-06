#ifndef TRACER_H_
#define TRACER_H_

#include "gloo/Scene.hpp"
#include "gloo/Material.hpp"
#include "gloo/lights/LightBase.hpp"
#include "gloo/components/LightComponent.hpp"

#include "Ray.hpp"
#include "HitRecord.hpp"
#include "TracingComponent.hpp"
#include "CubeMap.hpp"
#include "PerspectiveCamera.hpp"
#include "Camera.hpp"
#include "SpiralCamera.hpp"

namespace GLOO {

enum class FilterType {
  Box,   // Uniform weighting (equal weights for all samples)
  Tent   // Tent filter (higher weight near pixel center)
};

class Tracer {
 public:
  Tracer(Camera* camera,
         const glm::ivec2& image_size,
         size_t max_bounces,
         const glm::vec3& background_color,
         const CubeMap* cube_map,
         bool shadows_enabled)
      : camera_(camera),
        image_size_(image_size),
        max_bounces_(max_bounces),
        samples_per_dim_(1),  // Default to 1 sample per pixel (no AA - matches original behavior)
        filter_type_(FilterType::Box),  // Default to box filter (equal weights)
        background_color_(background_color),
        cube_map_(cube_map),
        shadows_enabled_(shadows_enabled),
        scene_ptr_(nullptr) {
  }
  void Render(const Scene& scene, const std::string& output_file);

  // Setters for anti-aliasing settings
  void SetSamplesPerDim(int samples_per_dim);
  void SetFilterType(FilterType filter_type);

 private:
  glm::vec3 TraceRay(const Ray& ray, size_t bounces, HitRecord& record) const;

  glm::vec3 GetBackgroundColor(const glm::vec3& direction) const;

  // Compute filter weight for a sample at relative offset (dx, dy) within pixel
  float ComputeFilterWeight(float dx, float dy, FilterType filter_type) const;

  Camera* camera_;  // Pointer to camera (owned by caller, or use unique_ptr)
  glm::ivec2 image_size_;
  size_t max_bounces_;
  int samples_per_dim_;  // Number of samples per pixel dimension (samples_per_dim_^2 total samples per pixel)
  FilterType filter_type_;  // Filter type for supersampling

  std::vector<TracingComponent*> tracing_components_;
  std::vector<LightComponent*> light_components_;
  glm::vec3 background_color_;
  const CubeMap* cube_map_;
  bool shadows_enabled_;

  const Scene* scene_ptr_;
};
}  // namespace GLOO

#endif
