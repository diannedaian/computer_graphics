#include "SkeletonNode.hpp"

#include <fstream>
#include <sstream>
#include "gloo/utils.hpp"
#include "gloo/InputManager.hpp"
#include "gloo/MeshLoader.hpp"
#include "gloo/debug/PrimitiveFactory.hpp"
#include "gloo/shaders/PhongShader.hpp"
#include "gloo/components/ShadingComponent.hpp"
#include "gloo/components/RenderingComponent.hpp"

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

// void SkeletonNode::DecorateTree() {
//   // TODO: set up addtional nodes, add necessary components here.
//   // You should create one set of nodes/components for skeleton mode
//   // (spheres for joints and cylinders for bones), and another set for
//   // SSD mode (you could just use a single node with a RenderingComponent
//   // that is linked to a VertexObject with the mesh information. Then you
//   // only need to update the VertexObject - updating vertex positions and
//   // recalculating the normals, etc.).

//   // The code snippet below shows how to add a sphere node to a joint.
//   // Suppose you have created member variables shader_ of type
//   // std::shared_ptr<PhongShader>, and sphere_mesh_ of type
//   // std::shared_ptr<VertexObject>.
//   // Here sphere_nodes_ptrs_ is a std::vector<SceneNode*> that stores the
//   // pointer so the sphere nodes can be accessed later to change their
//   // positions. joint_ptr is assumed to be one of the joint node you created
//   // from LoadSkeletonFile (e.g. you've stored a std::vector<SceneNode*> of
//   // joint nodes as a member variable and joint_ptr is one of the elements).
//   //
//   // auto sphere_node = make_unique<SceneNode>();
//   // sphere_node->CreateComponent<ShadingComponent>(shader_);
//   // sphere_node->CreateComponent<RenderingComponent>(sphere_mesh_);
//   // sphere_nodes_ptrs_.push_back(sphere_node.get());
//   // joint_ptr->AddChild(std::move(sphere_node));
//   shader_ = std::make_shared<PhongShader>();

//   // reusable mesh
//   float joint_radius = 0.03f;
//   sphere_mesh_ = PrimitiveFactory::CreateSphere(joint_radius, 15, 15);
//   float bone_radius = 0.0015f;
//   cylinder_mesh_ = PrimitiveFactory::CreateCylinder(bone_radius, 1.0f, 10);

//   for (size_t i = 0; i < joint_nodes_.size(); ++i) {
//     SceneNode* joint_ptr = joint_nodes_[i];
//     auto sphere_node = make_unique<SceneNode>();
//     sphere_node->CreateComponent<ShadingComponent>(shader_);
//     sphere_node->CreateComponent<RenderingComponent>(sphere_mesh_);
//     sphere_nodes_.push_back(sphere_node.get());
//     joint_ptr->AddChild(std::move(sphere_node));

//     if (joint_ptr->GetParentPtr() != this) {


//   }


// }
void SkeletonNode::DecorateTree() {
  // --- A. Setup Shared Resources ---
  shader_ = std::make_shared<PhongShader>();

  float joint_radius = 0.03f;
  sphere_mesh_ = PrimitiveFactory::CreateSphere(joint_radius, 15, 15);

  float bone_radius = 0.015f;
  // Base cylinder has radius 0.015 and height 1.0, positioned along +Y axis, bottom at (0,0,0).
  cylinder_mesh_ = PrimitiveFactory::CreateCylinder(bone_radius, 1.0f, 10);

  // --- B. Decorate Joints with Spheres and Bones ---

  // Iterate through all joints.
  for (size_t i = 0; i < joint_nodes_.size(); ++i) {
    SceneNode* joint_ptr = joint_nodes_[i];

    // 1. Add Sphere for the Joint (Always centered at the joint's origin)
    auto sphere_node = make_unique<SceneNode>();
    sphere_node->CreateComponent<ShadingComponent>(shader_);
    sphere_node->CreateComponent<RenderingComponent>(sphere_mesh_);
    sphere_nodes_.push_back(sphere_node.get());
    joint_ptr->AddChild(std::move(sphere_node));

    // 2. Add Cylinder (Bone) between this joint (Child) and its Parent
    // Condition to skip the skeleton root: The root's parent is 'this' (the SkeletonNode),
    // which is not a joint and should not have a bone extending to it.
    // The skeleton root joint is the first element, so we can also check i > 0.
    if (i > 0) {

      SceneNode* parent_ptr = joint_ptr->GetParentPtr();

      // Create the bone node and attach it as a child of the PARENT JOINT.
      auto bone_node = make_unique<SceneNode>();
      SceneNode* bone_node_ptr = bone_node.get();
      bone_node->CreateComponent<ShadingComponent>(shader_);
      bone_node->CreateComponent<RenderingComponent>(cylinder_mesh_);

      parent_ptr->AddChild(std::move(bone_node));

      cylinder_nodes_.push_back(bone_node_ptr);

      // --- Bone Transformation Logic (Relative to Parent) ---

      // Get the position of the child joint relative to the parent (its local position).
      glm::vec3 child_local_pos = joint_ptr->GetTransform().GetPosition();

      float length = glm::length(child_local_pos);

      if (length < 1e-6) continue;

      glm::vec3 direction = glm::normalize(child_local_pos);

      // Cylinder is built along +Y axis
      glm::vec3 y_axis(0.f, 1.f, 0.f);

      // Find the rotation required to align Y-axis with the 'direction' vector.
      glm::vec3 rotation_axis = glm::cross(y_axis, direction);
      float angle = acosf(glm::dot(y_axis, direction));

      Transform& bone_transform = bone_node_ptr->GetTransform();

      // 1. No Translation: The bone is already starting at (0, 0, 0), which is the parent's position.
      // bone_transform.SetPosition(glm::vec3(0.f)); // Implicitly zero/identity

      // 2. Rotation: Align the bone with the child joint.
      if (angle > 1e-6) {
          bone_transform.SetRotation(glm::normalize(rotation_axis), angle);
      }

      // 3. Scaling: Stretch the bone along Y to match the length.
      // The local bone height is 1.0, scaling by 'length' gives the correct bone length.
      bone_transform.SetScale(glm::vec3(1.f, length, 1.f));

      // Note: No translation is needed because the cylinder's bottom (0,0,0) is
      // where the parent joint is, and the cylinder's top (0, length, 0) is
      // where the child joint is (after rotation and scaling).
    }
  }
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
  std::ifstream file(path);
  std::string line;
  std::vector<SceneNode*> joint_ptrs;  // Store raw pointers for hierarchy building

  while (std::getline(file, line)) {
    std::istringstream iss(line);
    float x, y, z;           // Translation coordinates relative to parent
    int parent_index;        // Parent joint index (-1 for root)

    iss >> x >> y >> z >> parent_index;

    // Create a new joint node
    auto joint_node = make_unique<SceneNode>();
    joint_node->GetTransform().SetPosition(glm::vec3(x, y, z));

    // Store pointer for hierarchy building
    SceneNode* joint_ptr = joint_node.get();
    joint_ptrs.push_back(joint_ptr);
    joint_nodes_.push_back(joint_ptr);

    // Add to hierarchy
    if (parent_index == -1) {
      // Root joint - add as child of this SkeletonNode
      AddChild(std::move(joint_node));
    } else {
      // Child joint - add to its parent
      joint_ptrs[parent_index]->AddChild(std::move(joint_node));
    }
  }
  file.close();
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
