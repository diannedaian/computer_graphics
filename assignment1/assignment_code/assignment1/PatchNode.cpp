#include "PatchNode.hpp"

#include "gloo/components/RenderingComponent.hpp"
#include "gloo/components/ShadingComponent.hpp"
#include "gloo/components/MaterialComponent.hpp"
#include "gloo/shaders/PhongShader.hpp"
#include "gloo/Material.hpp"

namespace GLOO {
PatchNode::PatchNode(SplineBasis basis, const std::vector<glm::vec3>& control_points, const glm::vec3& color) {
  shader_ = std::make_shared<PhongShader>();
  patch_mesh_ = std::make_shared<VertexObject>();

  // TODO: this node should represent a single tensor product patch.
  // Think carefully about what data defines a patch and how you can
  // render it.
  spline_basis_ = basis;
  control_points_ = control_points;
  patch_color_ = color;

  InitPatch();
  PlotPatch();

}

void PatchNode::InitPatch() {
  auto patch_node = make_unique<SceneNode>();
  auto shader = std::make_shared<PhongShader>();
  patch_node->CreateComponent<ShadingComponent>(shader);

  // Create a material with the specified color
  auto material = std::make_shared<Material>(
    patch_color_ * 0.2f,          // ambient: darker version
    patch_color_,                 // diffuse: main color
    patch_color_ * 1.2f,          // specular: brighter version
    32.0f                        // shininess
  );
  patch_node->CreateComponent<MaterialComponent>(material);

  PlotPatch();
  auto& rc = patch_node->CreateComponent<RenderingComponent>(patch_mesh_);
  rc.SetDrawMode(DrawMode::Triangles);
  AddChild(std::move(patch_node));
}

PatchPoint PatchNode::EvalPatch(float u, float v) {
  // Clamp parameters to [0,1] range
  u = glm::clamp(u, 0.0f, 1.0f);
  v = glm::clamp(v, 0.0f, 1.0f);

  // Check if we have enough control points
  if (control_points_.size() < 16) {
    std::cout << "Warning: Not enough control points for patch evaluation. Have "
              << control_points_.size() << ", need 16." << std::endl;
    return PatchPoint{glm::vec3(0.0f), glm::vec3(0.0f, 0.0f, 1.0f), 0.0f};
  }

  // Create power basis vectors
  glm::vec4 U(1.0f, u, u*u, u*u*u);
  glm::vec4 V(1.0f, v, v*v, v*v*v);

  // Create derivative vectors
  glm::vec4 dU(0.0f, 1.0f, 2.0f*u, 3.0f*u*u);
  glm::vec4 dV(0.0f, 1.0f, 2.0f*v, 3.0f*v*v);

  // Create second derivative vectors
  glm::vec4 ddU(0.0f, 0.0f, 2.0f, 6.0f*u);
  glm::vec4 ddV(0.0f, 0.0f, 2.0f, 6.0f*v);

  // Select basis matrix
  glm::mat4 B;
  if (spline_basis_ == SplineBasis::Bezier) {
    B = glm::mat4(1, 0, 0, 0,
                  -3, 3, 0, 0,
                  3, -6, 3, 0,
                  -1, 3, -3, 1);
  } else {
    B = (1.0f / 6.0f) * glm::mat4(1, 4, 1, 0,
                                  -3, 0, 3, 0,
                                  3, -6, 3, 0,
                                  -1, 3, -3, 1);
  }

  // Compute blending weights
  glm::vec4 alpha = B * U;
  glm::vec4 beta = B * V;

  // Compute derivative weights
  glm::vec4 dalpha = B * dU;
  glm::vec4 dbeta = B * dV;

  // Compute second derivative weights
  glm::vec4 ddalpha = B * ddU;
  glm::vec4 ddbeta = B * ddV;

  // Initialize surface point and partial derivatives
  glm::vec3 P(0.0f);
  glm::vec3 Su(0.0f);
  glm::vec3 Sv(0.0f);
  glm::vec3 Suu(0.0f);
  glm::vec3 Svv(0.0f);
  glm::vec3 Suv(0.0f);

  // Evaluate surface point and partial derivatives using tensor product formula
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      int index = i * 4 + j;
      glm::vec3 control_point = control_points_[index];
      float weight = alpha[i] * beta[j];
      float dweight_u = dalpha[i] * beta[j];
      float dweight_v = alpha[i] * dbeta[j];
      float ddweight_uu = ddalpha[i] * beta[j];
      float ddweight_vv = alpha[i] * ddbeta[j];
      float ddweight_uv = dalpha[i] * dbeta[j];

      P += weight * control_point;
      Su += dweight_u * control_point;
      Sv += dweight_v * control_point;
      Suu += ddweight_uu * control_point;
      Svv += ddweight_vv * control_point;
      Suv += ddweight_uv * control_point;
    }
  }

  // Calculate surface normal
  glm::vec3 N = glm::normalize(glm::cross(Su, Sv));

  // calculate coefficients of the First Fundamental Form (I)
  float E = glm::dot(Su, Su);
  float F = glm::dot(Su, Sv);
  float G = glm::dot(Sv, Sv);

  // 2. Calculate coefficients of the Second Fundamental Form (II)

  float e = glm::dot(Suu, N);
  float f = glm::dot(Suv, N);
  float g = glm::dot(Svv, N);

  // 3. Calculate Gaussian curvature K = (eg - f^2) / (EG - F^2)
  float numerator = e * g - f * f;
  float denominator = E * G - F * F;

  float curvature;
  if (denominator < 1e-6) { // Check for singularity
    curvature = 0.0f;
  } else {
    curvature = numerator / denominator;
  }

  return PatchPoint{P, N, curvature};
}

// glm::vec3 PatchNode::CurvatureToColor(float curvature) {
//   // Find the range of curvature values and map them to a full color spectrum
//   float K_SCALE = 0.2f;
//   float scaled_curvature = curvature * K_SCALE;

//   float t = glm::clamp(scaled_curvature, 0.0f, 1.0f);


//   if (t < 0.33f) {
//     // Low to medium curvature: Blue to Green
//     float local_t = t / 0.33f; // Map [0, 0.33] to [0, 1]
//     return glm::vec3(0.0f, local_t, 1.0f - local_t); // Blue to Green
//   } else if (t < 0.66f) {
//     // Medium to high curvature: Green to Yellow
//     float local_t = (t - 0.33f) / 0.33f; // Map [0.33, 0.66] to [0, 1]
//     return glm::vec3(local_t, 1.0f, 0.0f); // Green to Yellow
//   } else {
//     // High curvature: Yellow to Red
//     float local_t = (t - 0.66f) / 0.34f; // Map [0.66, 1.0] to [0, 1]
//     return glm::vec3(1.0f, 1.0f - local_t, 0.0f); // Yellow to Red
//   }
// }

void PatchNode::PlotPatch() {
  std::cout << "PlotPatch: Starting with " << control_points_.size() << " control points" << std::endl;

  auto positions = make_unique<PositionArray>();
  auto normals = make_unique<NormalArray>();
  auto indices = make_unique<IndexArray>();

  // Sample the surface at regular intervals
  for (int i = 0; i <= N_SUBDIV_; i++) {
    for (int j = 0; j <= N_SUBDIV_; j++) {
      float u = static_cast<float>(i) / N_SUBDIV_;
      float v = static_cast<float>(j) / N_SUBDIV_;

      PatchPoint point = EvalPatch(u, v);
      positions->push_back(point.P);
      normals->push_back(point.N);
    }
  }

  std::cout << "PlotPatch: Generated " << positions->size() << " vertices" << std::endl;

  // Create triangular mesh
  for (int i = 0; i < N_SUBDIV_; i++) {
    for (int j = 0; j < N_SUBDIV_; j++) {
      int id1 = i * (N_SUBDIV_ + 1) + j;
      int id2 = id1 + 1;
      int id3 = (i + 1) * (N_SUBDIV_ + 1) + j;
      int id4 = id3 + 1;

      // Triangle 1 and 2
      indices->push_back(id1);
      indices->push_back(id2);
      indices->push_back(id3);
      indices->push_back(id2);
      indices->push_back(id3);
      indices->push_back(id4);
    }
  }

  patch_mesh_->UpdatePositions(std::move(positions));
  patch_mesh_->UpdateNormals(std::move(normals));
  patch_mesh_->UpdateIndices(std::move(indices));
}
}  // namespace GLOO
