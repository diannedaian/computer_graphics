#include "MeshViewerApp.hpp"
#include "TeapotNode.hpp"
#include "Lightmove.hpp"
#include "DirectionalLightToggle.hpp"

#include "glm/gtx/string_cast.hpp"

#include "gloo/components/CameraComponent.hpp"
#include "gloo/components/LightComponent.hpp"
#include "gloo/lights/PointLight.hpp"
#include "gloo/lights/AmbientLight.hpp"
#include "gloo/cameras/BasicCameraNode.hpp"
#include "gloo/lights/DirectionalLight.hpp"

namespace GLOO {
MeshViewerApp::MeshViewerApp(const std::string& app_name,
                             glm::ivec2 window_size)
    : Application(app_name, window_size) {
}

void MeshViewerApp::SetupScene() {
  SceneNode& root = scene_->GetRootNode();

  auto camera_node = make_unique<BasicCameraNode>();
  camera_node->GetTransform().SetPosition(glm::vec3(0.0f, 1.5f, 10.0f));
  scene_->ActivateCamera(camera_node->GetComponentPtr<CameraComponent>());
  root.AddChild(std::move(camera_node));

  auto ambient_light = std::make_shared<AmbientLight>();
  ambient_light->SetAmbientColor(glm::vec3(0.2f));
  root.CreateComponent<LightComponent>(ambient_light);

  auto point_light = std::make_shared<PointLight>();
  point_light->SetDiffuseColor(glm::vec3(0.8f, 0.8f, 0.8f));
  point_light->SetSpecularColor(glm::vec3(1.0f, 1.0f, 1.0f));
  point_light->SetAttenuation(glm::vec3(1.0f, 0.09f, 0.032f));
  auto point_light_node = make_unique<LightControlNode>();
  point_light_node->CreateComponent<LightComponent>(point_light);
  point_light_node->GetTransform().SetPosition(glm::vec3(0.0f, 4.0f, 5.f));
  root.AddChild(std::move(point_light_node));

  // Create directional light node
  auto directional_light = std::make_shared<DirectionalLight>();
  directional_light->SetDirection(glm::vec3(0.0f, -1.0f, -0.5f));
  directional_light->SetDiffuseColor(glm::vec3(0.8f, 0.8f, 0.8f));
  directional_light->SetSpecularColor(glm::vec3(1.0f, 1.0f, 1.0f));
  auto directional_light_node = make_unique<SceneNode>();
  directional_light_node->CreateComponent<LightComponent>(directional_light);

  // Create toggle controller for directional light
  auto directional_toggle = make_unique<DirectionalLightToggle>();
  directional_toggle->SetDirectionalLightNode(directional_light_node.get());

  root.AddChild(std::move(directional_light_node));
  root.AddChild(std::move(directional_toggle));

  // Add the teapot to the scene
  root.AddChild(make_unique<TeapotNode>());
}
}  // namespace GLOO
