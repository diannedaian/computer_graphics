#ifndef PATCH_NODE_H_
#define PATCH_NODE_H_

#include <string>
#include <vector>

#include "gloo/SceneNode.hpp"
#include "gloo/VertexObject.hpp"
#include "gloo/shaders/ShaderProgram.hpp"

#include "CurveNode.hpp"

namespace GLOO {
struct PatchPoint {
  glm::vec3 P;
  glm::vec3 N;
  float curvature;
};

class PatchNode : public SceneNode {
 public:
  PatchNode(SplineBasis basis, const std::vector<glm::vec3>& control_points, const glm::vec3& color = glm::vec3(0.0f, 0.0f, 0.8f));

 private:
  void PlotPatch();
  void InitPatch();
  PatchPoint EvalPatch(float u, float v);
  glm::vec3 CurvatureToColor(float curvature);

  std::vector<glm::mat4> Gs_;
  SplineBasis spline_basis_;
  std::vector<glm::vec3> control_points_;
  glm::vec3 patch_color_;

  std::shared_ptr<VertexObject> patch_mesh_;
  std::shared_ptr<ShaderProgram> shader_;

  const int N_SUBDIV_ = 50;
};
}  // namespace GLOO

#endif
