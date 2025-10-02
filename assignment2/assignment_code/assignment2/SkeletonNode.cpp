#include "SkeletonNode.hpp"

#include "gloo/utils.hpp"
#include "gloo/InputManager.hpp"
#include "gloo/MeshLoader.hpp"

namespace GLOO {
SkeletonNode::SkeletonNode(const std::string& filename)
    : SceneNode(), draw_mode_(DrawMode::Skeleton) {
  LoadAllFiles(filename);
  DecorateTree();

  // Force initial update.
  OnJointChanged();
}

void SkeletonNode::ToggleDrawMode() {
  draw_mode_ =
      draw_mode_ == DrawMode::Skeleton ? DrawMode::SSD : DrawMode::Skeleton;
  // TODO: implement here toggling between skeleton mode and SSD mode.
  // The current mode is draw_mode_;
  // Hint: you may find SceneNode::SetActive convenient here as
  // inactive nodes will not be picked up by the renderer.
}

void SkeletonNode::DecorateTree() {
  // TODO: set up addtional nodes, add necessary components here.
  // You should create one set of nodes/components for skeleton mode
  // (spheres for joints and cylinders for bones), and another set for
  // SSD mode (you could just use a single node with a RenderingComponent
  // that is linked to a VertexObject with the mesh information. Then you
  // only need to update the VertexObject - updating vertex positions and
  // recalculating the normals, etc.).

  // The code snippet below shows how to add a sphere node to a joint.
  // Suppose you have created member variables shader_ of type
  // std::shared_ptr<PhongShader>, and sphere_mesh_ of type
  // std::shared_ptr<VertexObject>.
  // Here sphere_nodes_ptrs_ is a std::vector<SceneNode*> that stores the
  // pointer so the sphere nodes can be accessed later to change their
  // positions. joint_ptr is assumed to be one of the joint node you created
  // from LoadSkeletonFile (e.g. you've stored a std::vector<SceneNode*> of
  // joint nodes as a member variable and joint_ptr is one of the elements).
  //
  // auto sphere_node = make_unique<SceneNode>();
  // sphere_node->CreateComponent<ShadingComponent>(shader_);
  // sphere_node->CreateComponent<RenderingComponent>(sphere_mesh_);
  // sphere_nodes_ptrs_.push_back(sphere_node.get());
  // joint_ptr->AddChild(std::move(sphere_node));
}

void SkeletonNode::Update(double delta_time) {
  // Prevent multiple toggle.
  static bool prev_released = true;
  if (InputManager::GetInstance().IsKeyPressed('S')) {
    if (prev_released) {
      ToggleDrawMode();
    }
    prev_released = false;
  } else if (InputManager::GetInstance().IsKeyReleased('S')) {
    prev_released = true;
  }
}

void SkeletonNode::OnJointChanged() {
  // TODO: this method is called whenever the values of UI sliders change.
  // The new Euler angles (represented as EulerAngle struct) can be retrieved
  // from linked_angles_ (a std::vector of EulerAngle*).
  // The indices of linked_angles_ align with the order of the joints in .skel
  // files. For instance, *linked_angles_[0] corresponds to the first line of
  // the .skel file.
}

void SkeletonNode::LinkRotationControl(const std::vector<EulerAngle*>& angles) {
  linked_angles_ = angles;
}

void SkeletonNode::LoadSkeletonFile(const std::string& path) {
  // TODO: load skeleton file and build the tree of joints.
}

void SkeletonNode::LoadMeshFile(const std::string& filename) {
  std::shared_ptr<VertexObject> vtx_obj =
      MeshLoader::Import(filename).vertex_obj;
  // TODO: store the bind pose mesh in your preferred way.
}

void SkeletonNode::LoadAttachmentWeights(const std::string& path) {
  // TODO: load attachment weights.
}

void SkeletonNode::LoadAllFiles(const std::string& prefix) {
  std::string prefix_full = GetAssetDir() + prefix;
  LoadSkeletonFile(prefix_full + ".skel");
  LoadMeshFile(prefix + ".obj");
  LoadAttachmentWeights(prefix_full + ".attach");
}
}  // namespace GLOO
