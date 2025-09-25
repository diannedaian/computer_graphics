#include "CurveNode.hpp"

#include "gloo/debug/PrimitiveFactory.hpp"
#include "gloo/components/RenderingComponent.hpp"
#include "gloo/components/ShadingComponent.hpp"
#include "gloo/components/MaterialComponent.hpp"
#include "gloo/shaders/PhongShader.hpp"
#include "gloo/shaders/SimpleShader.hpp"
#include "gloo/InputManager.hpp"

namespace GLOO {
CurveNode::CurveNode(SplineBasis basis, const std::vector<glm::vec3>& control_points, bool allow_toggle) {
  // TODO: this node should represent a single spline curve.
  // Think carefully about what data defines a curve and how you can
  // render it.

  // Initialize the VertexObjects and shaders used to render the control points,
  // the curve, and the tangent line.
  std::cout << "Curve node constructor called" << std::endl;
  allow_toggle_ = allow_toggle;
  sphere_mesh_ = PrimitiveFactory::CreateSphere(0.015f, 25, 25);
  curve_polyline_ = std::make_shared<VertexObject>();
  tangent_line_ = std::make_shared<VertexObject>();
  shader_ = std::make_shared<PhongShader>();
  polyline_shader_ = std::make_shared<SimpleShader>();
  spline_basis_ = basis;
  control_points_ = control_points;


  InitCurve();
  PlotCurve();
}

void CurveNode::Update(double delta_time) {
  // Prevent multiple toggles from a single key press.
  if (!allow_toggle_) {
    return;
  }
  static bool prev_released = true;

  if (InputManager::GetInstance().IsKeyPressed('T')) {
    if (prev_released) {
      // Toggle the evaluation basis (this also changes control point color)
      ToggleSplineBasis();
    }
    prev_released = false;
  } else if (InputManager::GetInstance().IsKeyPressed('B')) {
    if (prev_released) {
      // Convert geometry from Bezier to B-Spline.
      ConvertGeometry(true);
    }
    prev_released = false;
  } else if (InputManager::GetInstance().IsKeyPressed('Z')) {
    if (prev_released) {
      // Convert geometry from B-Spline to Bezier.
      ConvertGeometry(false);
    }
    prev_released = false;
  } else {
    prev_released = true;
  }
  // PlotCurve();
  // PlotTangentLine();
}

void CurveNode::ToggleSplineBasis() {
  // TODO: implement toggling between Bezier and B-Spline bases.
  if (spline_basis_ == SplineBasis::Bezier) {
    spline_basis_ = SplineBasis::BSpline;
    //change color of Bspline control point to green
    for (auto& control_node : control_nodes_) {
      auto material_comp = control_node->GetComponentPtr<MaterialComponent>();
      if (material_comp != nullptr) {
        auto& material = material_comp->GetMaterial();
        material.SetDiffuseColor(glm::vec3(0.0f, 1.0f, 0.0f));
        material.SetAmbientColor(glm::vec3(0.0f, 1.0f, 0.0f));
        material.SetSpecularColor(glm::vec3(0.0f, 1.0f, 0.0f));
      }
    }
  } else {
    spline_basis_ = SplineBasis::Bezier;
    //change color of Bezier control point to red
    for (auto& control_node : control_nodes_) {
      auto material_comp = control_node->GetComponentPtr<MaterialComponent>();
      if (material_comp != nullptr) {
        auto& material = material_comp->GetMaterial();
        material.SetDiffuseColor(glm::vec3(1.0f, 0.0f, 0.0f));
        material.SetAmbientColor(glm::vec3(1.0f, 0.0f, 0.0f));
        material.SetSpecularColor(glm::vec3(1.0f, 0.0f, 0.0f));
      }
    }
  }
  PlotCurve();
  PlotTangentLine();
}

void CurveNode::ConvertGeometry(bool bezier_to_bspline) {
  glm::mat4 Bsrc, Bdst;

  // Determine the source and destination bases based on the boolean parameter.
  if (bezier_to_bspline) {
    Bsrc = B_bezier;
    Bdst = B_bspline;
    std::cout << "Converting Bezier to B-Spline" << std::endl;
  } else {
    Bsrc = B_bspline;
    Bdst = B_bezier;
    std::cout << "Converting B-Spline to Bezier" << std::endl;
  }

  // The conversion matrix M is computed as the transpose of (B_dst * inverse(B_src)).
  glm::mat4 M = glm::transpose(Bdst * glm::inverse(Bsrc));

  // Gather control points into vectors for matrix multiplication.
  glm::vec4 X, Y, Z;
  for (int i = 0; i < 4; ++i) {
    const glm::vec3 p = control_nodes_[i]->GetTransform().GetPosition();
    X[i] = p.x; Y[i] = p.y; Z[i] = p.z;
  }

  // Apply the conversion matrix to each coordinate vector.
  X = M * X; Y = M * Y; Z = M * Z;

  // Update the positions of the control point nodes with the new geometry.
  for (int i = 0; i < 4; ++i) {
    control_nodes_[i]->GetTransform().SetPosition(glm::vec3(X[i], Y[i], Z[i]));
  }

  //spline_basis_ = dst;
  PlotCurve();
  PlotTangentLine();
}

CurvePoint CurveNode::EvalCurve(float t) {
  if (control_nodes_.empty()) {
    std::cout << "EvalCurve: No control nodes!" << std::endl;
    return CurvePoint();
  }

  t = glm::clamp(t, 0.0f, 1.0f);

  // Use the first 4 control nodes for evaluation
  glm::vec3 p0 = control_nodes_[0]->GetTransform().GetPosition();
  glm::vec3 p1 = control_nodes_[1]->GetTransform().GetPosition();
  glm::vec3 p2 = control_nodes_[2]->GetTransform().GetPosition();
  glm::vec3 p3 = control_nodes_[3]->GetTransform().GetPosition();

  // Set up the parameter vector T and its derivative dT
  glm::vec4 T(1, t, t*t, t*t*t);
  glm::vec4 dT(0, 1, 2*t, 3*t*t);

  glm::mat4 B;
  if (spline_basis_ == SplineBasis::Bezier) {
    B = B_bezier;
  } else {
    B = B_bspline;
  }

  // Calculate coefficients
  glm::vec4 coeffs = B * T;
  glm::vec4 deriv_coeffs = B * dT;

  // Linear combination of control points
  glm::vec3 point(0.0f);
  glm::vec3 tangent(0.0f);

  point += coeffs[0] * p0;
  point += coeffs[1] * p1;
  point += coeffs[2] * p2;
  point += coeffs[3] * p3;

  tangent += deriv_coeffs[0] * p0;
  tangent += deriv_coeffs[1] * p1;
  tangent += deriv_coeffs[2] * p2;
  tangent += deriv_coeffs[3] * p3;

  return CurvePoint{point, glm::normalize(tangent)};
}

void CurveNode::InitCurve() {
  // First, create the control point nodes.
  for (auto& pos : control_points_) {
    auto control_point_node = make_unique<SceneNode>();
    control_point_node->CreateComponent<ShadingComponent>(shader_);
    control_point_node->CreateComponent<RenderingComponent>(sphere_mesh_);

    glm::vec3 node_color;
    if (spline_basis_ == SplineBasis::Bezier) {
      node_color = glm::vec3(1.f, 0.f, 0.f); // Red for Bézier
    } else {
      node_color = glm::vec3(0.f, 1.f, 0.f); // Green for B-Spline
    }

    auto material = std::make_shared<Material>(
        node_color,
        node_color,
        glm::vec3(1.f, 0.f, 0.f),
        0.0f
    );
    control_point_node->CreateComponent<MaterialComponent>(material);
    control_point_node->GetTransform().SetPosition(pos);

    control_nodes_.push_back(control_point_node.get());
    AddChild(std::move(control_point_node));
  }

  // Create the curve line node
  curve_node_ = make_unique<SceneNode>();
  curve_node_->CreateComponent<ShadingComponent>(polyline_shader_);

  // Plot the curve first to populate the VertexObject
  PlotCurve();

  auto& curve_rc = curve_node_->CreateComponent<RenderingComponent>(curve_polyline_);
  curve_rc.SetDrawMode(DrawMode::Lines);

  glm::vec3 curve_color(1.f, 1.f, 0.f); // Yellow for the curve
  auto curve_material = std::make_shared<Material>(curve_color, curve_color, curve_color, 0);
  curve_node_->CreateComponent<MaterialComponent>(curve_material);
  AddChild(std::move(curve_node_));

  // Create the tangent line node
  tangent_line_node_ = make_unique<SceneNode>();
  tangent_line_node_->CreateComponent<ShadingComponent>(polyline_shader_);

  // Plot the tangent line first to populate the VertexObject
  PlotTangentLine();

  auto& tangent_rc = tangent_line_node_->CreateComponent<RenderingComponent>(tangent_line_);
  tangent_rc.SetDrawMode(DrawMode::Lines);

  glm::vec3 tangent_color(1.f, 1.f, 1.f); // White for the tangent line
  auto tangent_material = std::make_shared<Material>(tangent_color, tangent_color, tangent_color, 0);
  tangent_line_node_->CreateComponent<MaterialComponent>(tangent_material);
  AddChild(std::move(tangent_line_node_));
}

void CurveNode::PlotCurve() {
  std::cout << "PlotCurve called" << std::endl;
  if (control_nodes_.empty()) {
    std::cout << "PlotCurve: No control nodes!" << std::endl;
    return;
  }

  auto positions = make_unique<PositionArray>();
  auto indices = make_unique<IndexArray>();

  // Generate curve points
  for (int i = 0; i <= N_SUBDIV_; i++) {
    float t = static_cast<float>(i) / N_SUBDIV_;
    CurvePoint point = EvalCurve(t);
    positions->push_back(point.P);
  }

  // Create line segments
  for (size_t i = 0; i < positions->size() - 1; i++) {
    indices->push_back(i);
    indices->push_back(i + 1);
  }

  // Update the VertexObject
  curve_polyline_->UpdatePositions(std::move(positions));
  curve_polyline_->UpdateIndices(std::move(indices));
}

void CurveNode::PlotControlPoints() {

  control_nodes_.clear();

  // Create and plot new control point nodes.
  for (size_t i = 0; i < control_points_.size(); i++) {
    auto point_node = make_unique<SceneNode>();
    point_node->GetTransform().SetPosition(control_points_[i]);

    auto& rc = point_node->CreateComponent<RenderingComponent>(sphere_mesh_);
    rc.SetDrawMode(DrawMode::Triangles);

    point_node->CreateComponent<ShadingComponent>(shader_);

    // Set the color based on the current spline basis.
    glm::vec3 color = (spline_basis_ == SplineBasis::Bezier) ?
      glm::vec3(1.0f, 0.0f, 0.0f) : glm::vec3(0.0f, 1.0f, 0.0f);

    auto material = std::make_shared<Material>(color, color, color, 0);
    point_node->CreateComponent<MaterialComponent>(material);

    // Store a raw pointer for later access and add the node to the scene.
    control_nodes_.push_back(point_node.get());
    AddChild(std::move(point_node));
  }
}

void CurveNode::PlotTangentLine() {
  if (control_nodes_.empty()) {
    return;
  }

  auto positions = make_unique<PositionArray>();
  auto indices = make_unique<IndexArray>();

  float t = 0.5f;
  CurvePoint point = EvalCurve(t);
  glm::vec3 p = point.P;
  glm::vec3 T = point.T;

  positions->push_back(p - 0.1f * T); // Start
  positions->push_back(p + 0.1f * T); // End 

  indices->push_back(0);
  indices->push_back(1);

  tangent_line_->UpdatePositions(std::move(positions));
  tangent_line_->UpdateIndices(std::move(indices));
}
}  // namespace GLOO
