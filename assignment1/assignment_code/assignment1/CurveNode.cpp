#include "CurveNode.hpp"

#include "gloo/debug/PrimitiveFactory.hpp"
#include "gloo/components/RenderingComponent.hpp"
#include "gloo/components/ShadingComponent.hpp"
#include "gloo/components/MaterialComponent.hpp"
#include "gloo/shaders/PhongShader.hpp"
#include "gloo/shaders/SimpleShader.hpp"
#include "gloo/InputManager.hpp"

namespace GLOO {
CurveNode::CurveNode() {
  // TODO: this node should represent a single spline curve.
  // Think carefully about what data defines a curve and how you can
  // render it.

  // Initialize the VertexObjects and shaders used to render the control points,
  // the curve, and the tangent line.
  sphere_mesh_ = PrimitiveFactory::CreateSphere(0.015f, 25, 25);
  curve_polyline_ = std::make_shared<VertexObject>();
  tangent_line_ = std::make_shared<VertexObject>();
  shader_ = std::make_shared<PhongShader>();
  polyline_shader_ = std::make_shared<SimpleShader>();

  InitCurve();
  PlotCurve();
}

void CurveNode::Update(double delta_time) {

  // Prevent multiple toggle.
  static bool prev_released = true;

  if (InputManager::GetInstance().IsKeyPressed('T')) {
    if (prev_released) {
      // TODO: implement toggling spline bases.
    }
    prev_released = false;
  } else if (InputManager::GetInstance().IsKeyPressed('B')) {
    if (prev_released) {
      // TODO: implement converting conrol point geometry from Bezier to
      // B-Spline basis.
    }
    prev_released = false;
  } else if (InputManager::GetInstance().IsKeyPressed('Z')) {
    if (prev_released) {
      // TODO: implement converting conrol point geometry from B-Spline to
      // Bezier basis.
    }
    prev_released = false;
  } else {
    prev_released = true;
  }
}

void CurveNode::ToggleSplineBasis() {
  // TODO: implement toggling between Bezier and B-Spline bases.
}

void CurveNode::ConvertGeometry() {
  // TODO: implement converting the control points between bases.
}

CurvePoint CurveNode::EvalCurve(float t) {
  // TODO: implement evaluating the spline curve at parameter value t.
  return CurvePoint();
}

void CurveNode::InitCurve() {
  // TODO: create all of the  nodes and components necessary for rendering the
  // curve, its control points, and its tangent line. You will want to use the
  // VertexObjects and shaders that are initialized in the class constructor.
}

void CurveNode::PlotCurve() {
  // TODO: plot the curve by updating the positions of its VertexObject.
}

void CurveNode::PlotControlPoints() {
  // TODO: plot the curve control points.
}

void CurveNode::PlotTangentLine() {
  // TODO: implement plotting a line tangent to the curve.
  // Below is a sample implementation for rendering a line segment
  // onto the screen. Note that this is just an example. This code
  // currently has nothing to do with the spline.

  auto line = std::make_shared<VertexObject>();

  auto positions = make_unique<PositionArray>();
  positions->push_back(glm::vec3(1.f, 0, -1.f));
  positions->push_back(glm::vec3(5.f, 5.f, 5.f));

  auto indices = make_unique<IndexArray>();
  indices->push_back(0);
  indices->push_back(1);

  line->UpdatePositions(std::move(positions));
  line->UpdateIndices(std::move(indices));

  auto shader = std::make_shared<SimpleShader>();

  auto line_node = make_unique<SceneNode>();
  line_node->CreateComponent<ShadingComponent>(shader);

  auto& rc = line_node->CreateComponent<RenderingComponent>(line);
  rc.SetDrawMode(DrawMode::Lines);

  glm::vec3 color(1.f, 1.f, 1.f);
  auto material = std::make_shared<Material>(color, color, color, 0);
  line_node->CreateComponent<MaterialComponent>(material);

  AddChild(std::move(line_node));
}
}  // namespace GLOO
