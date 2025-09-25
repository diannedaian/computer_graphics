#ifndef CURVE_NODE_H_
#define CURVE_NODE_H_

#include <string>
#include <vector>

#include "gloo/SceneNode.hpp"
#include "gloo/VertexObject.hpp"
#include "gloo/shaders/ShaderProgram.hpp"

namespace GLOO {

enum class SplineBasis { Bezier, BSpline };

struct CurvePoint {
  glm::vec3 P;
  glm::vec3 T;
};

class CurveNode : public SceneNode {
 public:
  CurveNode(SplineBasis basis, const std::vector<glm::vec3>& control_points, bool allow_toggle);
  void Update(double delta_time) override;

  std::vector<glm::vec3> control_points_;
  bool allow_toggle_;
  SplineBasis spline_basis_;
  void InitCurve();

 private:
  void ToggleSplineBasis();
  void ConvertGeometry(bool bezier_to_bspline);
  CurvePoint EvalCurve(float t);
  void PlotCurve();
  void PlotControlPoints();
  void PlotTangentLine();

  std::vector<SceneNode*> control_nodes_;
  std::unique_ptr<SceneNode> curve_node_;
  std::unique_ptr<SceneNode> tangent_line_node_;

  std::shared_ptr<VertexObject> sphere_mesh_;
  std::shared_ptr<VertexObject> curve_polyline_;
  std::shared_ptr<VertexObject> tangent_line_;

  std::shared_ptr<ShaderProgram> shader_;
  std::shared_ptr<ShaderProgram> polyline_shader_;


  const int N_SUBDIV_ = 50;

  glm::mat4 B_bezier = glm::mat4(1,0,0,0,
    -3,3,0,0,
    3,-6,3,0,
    -1,3,-3,1
  );

  glm::mat4 B_bspline = (1.0f / 6.0f) * glm::mat4(1,4,1,0,
    -3,0,3,0,
    3,-6,3,0,
    -1,3,-3,1
  );
};
}  // namespace GLOO

#endif
