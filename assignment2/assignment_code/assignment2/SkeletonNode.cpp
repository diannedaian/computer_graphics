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
#include "gloo/alias_types.hpp"

namespace GLOO {
SkeletonNode::SkeletonNode(const std::string& filename)
    : SceneNode(), active_char_idx_(0), draw_mode_(DrawMode::Skeleton) {
  // Load shared resources
  LoadAllFiles(filename);

  // Initialize two characters at different positions
  InitializeCharacter(filename, glm::vec3(-1.0f, 0.0f, 0.0f));  // Character 1
  InitializeCharacter(filename, glm::vec3(1.0f, 0.0f, 0.0f));   // Character 2

  // Force initial update.
  OnJointChanged();
}

void SkeletonNode::ToggleDrawMode() {
  draw_mode_ =
      draw_mode_ == DrawMode::Skeleton ? DrawMode::SSD : DrawMode::Skeleton;

  // Apply mode change to ALL characters
  for (auto& character : characters_) {
    if (draw_mode_ == DrawMode::Skeleton) {
      // Show skeleton mode: activate sphere and cylinder nodes
      for (auto* sphere_node : character.sphere_nodes) {
        sphere_node->SetActive(true);
      }
      for (auto* cylinder_node : character.cylinder_nodes) {
        cylinder_node->SetActive(true);
      }
      // Hide SSD mesh
      if (character.ssd_node) {
        character.ssd_node->SetActive(false);
      }
    } else {
      // Show SSD mode: hide skeleton, show deformed mesh
      for (auto* sphere_node : character.sphere_nodes) {
        sphere_node->SetActive(false);
      }
      for (auto* cylinder_node : character.cylinder_nodes) {
        cylinder_node->SetActive(false);
      }
      // Show SSD mesh
      if (character.ssd_node) {
        character.ssd_node->SetActive(true);
      }
    }
  }
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
void SkeletonNode::DecorateTree(CharacterInstance& character) {
  // Iterate through all joints.
  for (size_t i = 0; i < character.joint_nodes.size(); ++i) {
    SceneNode* joint_ptr = character.joint_nodes[i];

    // 1. Add Sphere for the Joint (Always centered at the joint's origin)
    auto sphere_node = make_unique<SceneNode>();
    sphere_node->CreateComponent<ShadingComponent>(shader_);
    sphere_node->CreateComponent<RenderingComponent>(sphere_mesh_);
    character.sphere_nodes.push_back(sphere_node.get());
    joint_ptr->AddChild(std::move(sphere_node));

    // 2. Add Cylinder (Bone) between this joint (Child) and its Parent
    // Only create cylinders between actual joint nodes, not between joint and character root
    // Skip root joint (i = 0) and any joint whose parent is character root
    if (i > 0 && joint_ptr->GetParentPtr() != character.root_node.get()) {
      SceneNode* parent_ptr = joint_ptr->GetParentPtr();

      // Create the bone node and attach it as a child of the PARENT JOINT.
      auto bone_node = make_unique<SceneNode>();
      SceneNode* bone_node_ptr = bone_node.get();
      bone_node->CreateComponent<ShadingComponent>(shader_);
      bone_node->CreateComponent<RenderingComponent>(cylinder_mesh_);

      parent_ptr->AddChild(std::move(bone_node));

      character.cylinder_nodes.push_back(bone_node_ptr);

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

      // 2. Rotation: Align the bone with the child joint.
      if (angle > 1e-6) {
          bone_transform.SetRotation(glm::normalize(rotation_axis), angle);
      }

      // 3. Scaling: Stretch the bone along Y to match the length.
      bone_transform.SetScale(glm::vec3(1.f, length, 1.f));
    }
  }

  // --- C. Setup SSD Node ---
  // Create a node for the deformed mesh (initially inactive)
  auto ssd_node = make_unique<SceneNode>();
  character.ssd_node = ssd_node.get();
  character.ssd_node->CreateComponent<ShadingComponent>(shader_);
  // The RenderingComponent will be set up when the deformed mesh is created
  character.root_node->AddChild(std::move(ssd_node));

  // Initially hide SSD mode (show skeleton by default)
  character.ssd_node->SetActive(false);
}

void SkeletonNode::Update(double delta_time) {
  // Prevent multiple toggles with static flags
  static bool s_prev_released = true;
  static bool one_prev_released = true;
  static bool two_prev_released = true;

  // Toggle draw mode with 'S' key
  if (InputManager::GetInstance().IsKeyPressed('S')) {
    if (s_prev_released) {
      ToggleDrawMode();
    }
    s_prev_released = false;
  } else if (InputManager::GetInstance().IsKeyReleased('S')) {
    s_prev_released = true;
  }

  // Switch to character 1 with '1' key
  if (InputManager::GetInstance().IsKeyPressed('1')) {
    if (one_prev_released && active_char_idx_ != 0) {
      UpdateActiveCharacterPose(0);
    }
    one_prev_released = false;
  } else if (InputManager::GetInstance().IsKeyReleased('1')) {
    one_prev_released = true;
  }

  // Switch to character 2 with '2' key
  if (InputManager::GetInstance().IsKeyPressed('2')) {
    if (two_prev_released && active_char_idx_ != 1) {
      UpdateActiveCharacterPose(1);
    }
    two_prev_released = false;
  } else if (InputManager::GetInstance().IsKeyReleased('2')) {
    two_prev_released = true;
  }
}

// void SkeletonNode::OnJointChanged() {
//   // TODO: this method is called whenever the values of UI sliders change.
//   // The new Euler angles (represented as EulerAngle struct) can be retrieved
//   // from linked_angles_ (a std::vector of EulerAngle*).
//   // The indices of linked_angles_ align with the order of the joints in .skel
//   // files. For instance, *linked_angles_[0] corresponds to the first line of
//   // the .skel file.
//     // Setup shared resources (shader, sphere mesh for joints, cylinder mesh for bones).
//     shader_ = std::make_shared<PhongShader>();

//     float joint_radius = 0.03f;
//     sphere_mesh_ = PrimitiveFactory::CreateSphere(joint_radius, 15, 15);

//     float bone_radius = 0.015f;
//     cylinder_mesh_ = PrimitiveFactory::CreateCylinder(bone_radius, 1.0f, 10);

//     // Iterate over joints to add sphere meshes and bone cylinders.
//     for (size_t i = 0; i < joint_nodes_.size(); ++i) {
//       SceneNode* joint_ptr = joint_nodes_[i];

//       // Add sphere to represent the joint location.
//       auto sphere_node = make_unique<SceneNode>();
//       sphere_node->CreateComponent<ShadingComponent>(shader_);
//       sphere_node->CreateComponent<RenderingComponent>(sphere_mesh_);
//       sphere_nodes_.push_back(sphere_node.get());
//       joint_ptr->AddChild(std::move(sphere_node));

//       // If not the root joint, create a cylinder (bone) and attach it to the parent joint.
//       if (i > 0 && joint_ptr->GetParentPtr() != this) {
//         SceneNode* parent_ptr = joint_ptr->GetParentPtr();

//         auto bone_node = make_unique<SceneNode>();
//         SceneNode* bone_node_ptr = bone_node.get();
//         bone_node->CreateComponent<ShadingComponent>(shader_);
//         bone_node->CreateComponent<RenderingComponent>(cylinder_mesh_);
//         parent_ptr->AddChild(std::move(bone_node));
//         cylinder_nodes_.push_back(bone_node_ptr);

//         // Calculate initial bone length and rotation based on child's local position.
//         glm::vec3 child_local_pos = joint_ptr->GetTransform().GetPosition();
//         float length = glm::length(child_local_pos);

//         if (length < 1e-6) continue;

//         glm::vec3 direction = glm::normalize(child_local_pos);
//         glm::vec3 y_axis(0.f, 1.f, 0.f);

//         glm::vec3 rotation_axis = glm::cross(y_axis, direction);
//         float angle = acosf(glm::dot(y_axis, direction));

//         Transform& bone_transform = bone_node_ptr->GetTransform();

//         // Apply initial rotation and scale for the bone.
//         if (angle > 1e-6) {
//             bone_transform.SetRotation(glm::normalize(rotation_axis), angle);
//         }
//         bone_transform.SetScale(glm::vec3(1.f, length, 1.f));
//       }
//     }
//   }


  void SkeletonNode::OnJointChanged() {
    // Get the active character
    if (active_char_idx_ < 0 || active_char_idx_ >= static_cast<int>(characters_.size())) {
      return;
    }

    CharacterInstance& character = characters_[active_char_idx_];

    // Update joint rotations based on linked UI sliders.
    if (linked_angles_.size() != character.joint_nodes.size()) {
        return;
    }

    // Apply new Euler angle rotations to joints.
    for (size_t i = 0; i < linked_angles_.size(); ++i) {
      const EulerAngle& angle = *linked_angles_[i];
      SceneNode* joint_ptr = character.joint_nodes[i];

      // Use angles directly
      glm::vec3 euler_rads(angle.rx, angle.ry, angle.rz);

      // Convert to quaternion
      glm::quat rotation_quat(euler_rads);
      joint_ptr->GetTransform().SetRotation(rotation_quat);
    }

    // Recompute the scale and rotation of the bone cylinders to span the new joint positions.
    for (size_t i = 0; i < character.cylinder_nodes.size(); i++) {
      size_t child_joint_idx = i + 1;
      SceneNode* child_joint = character.joint_nodes[child_joint_idx];
      SceneNode* parent_joint = child_joint->GetParentPtr();

      if (parent_joint == character.root_node.get()) continue;

      // Use the child joint's local position relative to the parent to define the bone.
      glm::vec3 bone_dir_local = child_joint->GetTransform().GetPosition();
      float bone_length = glm::length(bone_dir_local);

      SceneNode* cylinder_node = character.cylinder_nodes[i];

      if (bone_length < 1e-6) {
        cylinder_node->SetActive(false);
        continue;
      }

      cylinder_node->SetActive(true);

      // Calculate rotation to align cylinder (along Y-axis) with the local bone direction.
      glm::vec3 cylinder_axis(0.0f, 1.0f, 0.0f);
      glm::vec3 direction = glm::normalize(bone_dir_local);
      glm::vec3 rotation_axis = glm::cross(cylinder_axis, direction);
      float rotation_angle = std::acos(glm::clamp(glm::dot(cylinder_axis, direction), -1.0f, 1.0f));

      Transform& transform = cylinder_node->GetTransform();

      // Set bone position (relative to parent joint) and length.
      transform.SetPosition(glm::vec3(0.0f, 0.0f, 0.0f));
      transform.SetScale(glm::vec3(1.0f, bone_length, 1.0f));

      // Apply the alignment rotation.
      if (glm::length(rotation_axis) > 1e-6) {
        rotation_axis = glm::normalize(rotation_axis);
        transform.SetRotation(rotation_axis, rotation_angle);
      } else {
        if (glm::dot(cylinder_axis, direction) < 0) {
          transform.SetRotation(glm::vec3(1.0f, 0.0f, 0.0f), glm::pi<float>());
        } else {
          transform.SetRotation(glm::quat(1.0f, 0.0f, 0.0f, 0.0f));
        }
      }
    }

    // Deform mesh for SSD mode
    DeformMesh(character);
}

void SkeletonNode::LinkRotationControl(const std::vector<EulerAngle*>& angles) {
  linked_angles_ = angles;
}

void SkeletonNode::LoadSkeletonFile(const std::string& path, CharacterInstance& character) {
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
    character.joint_nodes.push_back(joint_ptr);

    // Add to hierarchy
    if (parent_index == -1) {
      // Root joint - add as child of character's root node
      character.root_node->AddChild(std::move(joint_node));
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
  // Store the bind pose mesh for SSD mode
  bind_pose_mesh_ = vtx_obj;
}
void SkeletonNode::ComputeBindPoseMatrices(CharacterInstance& character) {
  // Compute bind pose transformation matrices Bi for each joint
  // Bi = joint's world transformation matrix in bind pose (LocalToWorld)
  character.bind_pose_matrices.clear();
  character.bind_pose_matrices.resize(character.joint_nodes.size());

  for (size_t i = 0; i < character.joint_nodes.size(); ++i) {
    SceneNode* joint = character.joint_nodes[i];

    // Use the correctly implemented GetLocalToWorldMatrix()
    character.bind_pose_matrices[i] = joint->GetTransform().GetLocalToWorldMatrix();
  }
}
// Inside GLOO::SkeletonNode
void SkeletonNode::ComputeNormals(const std::vector<glm::vec3>& vertices, std::vector<glm::vec3>& normals) {
  // Initialize all normals to zero
  std::fill(normals.begin(), normals.end(), glm::vec3(0.0f));

  // Get the indices from the original mesh
  const auto& indices = bind_pose_mesh_->GetIndices();

  // Process each triangle
  for (size_t i = 0; i < indices.size(); i += 3) {
    if (i + 2 >= indices.size()) break;

    unsigned int i0 = indices[i];
    unsigned int i1 = indices[i + 1];
    unsigned int i2 = indices[i + 2];

    if (i0 >= vertices.size() || i1 >= vertices.size() || i2 >= vertices.size()) {
      continue;
    }

    glm::vec3 v0 = vertices[i0];
    glm::vec3 v1 = vertices[i1];
    glm::vec3 v2 = vertices[i2];

    // Compute face normal (cross product of edges)
    glm::vec3 edge1 = v1 - v0;
    glm::vec3 edge2 = v2 - v0;
    glm::vec3 face_normal = glm::cross(edge1, edge2);

    // Compute face area (half the magnitude of the cross product is related to the magnitude of the cross product)
    float face_area = glm::length(face_normal); // Use the magnitude directly as the weight

    // Skip degenerate triangles
    if (face_area < 1e-8f) {
      continue;
    }

    glm::vec3 normalized_face_normal = face_normal / face_area;

    // Add weighted face normal to each vertex of the triangle
    // Weighting by the magnitude (2 * area) ensures area-proportional contribution
    normals[i0] += face_area * normalized_face_normal;
    normals[i1] += face_area * normalized_face_normal;
    normals[i2] += face_area * normalized_face_normal;
  }

  // Normalize all vertex normals
  for (size_t i = 0; i < normals.size(); ++i) {
    if (glm::length(normals[i]) > 1e-8f) {
      normals[i] = glm::normalize(normals[i]);
    } else {
      // Use default normal (e.g., up) if vector is zero
      normals[i] = glm::vec3(0.0f, 1.0f, 0.0f);
    }
  }
}


void SkeletonNode::DeformMesh(CharacterInstance& character) {
  if (!bind_pose_mesh_ || attachment_weights_.empty()) {
    return;
  }

  const auto& original_vertices = bind_pose_mesh_->GetPositions();
  std::vector<glm::vec3> deformed_vertices(original_vertices.size());

  // For each vertex, compute its deformed position
  for (size_t vertex_idx = 0; vertex_idx < original_vertices.size(); ++vertex_idx) {
    glm::vec3 original_pos = original_vertices[vertex_idx];
    glm::vec4 deformed_pos(0.0f);

    float total_weight = 0.0f;

    // 1. Accumulate transformations for joints 1 to m-1 (non-root joints)
    for (size_t joint_idx = 1; joint_idx < character.joint_nodes.size(); ++joint_idx) {
      // Attachment weights for vertex i are for joints 1, 2, ...
      float weight = 0.0f;
      if (vertex_idx < attachment_weights_.size() && joint_idx - 1 < attachment_weights_[vertex_idx].size()) {
        weight = attachment_weights_[vertex_idx][joint_idx - 1];
      }

      if (weight > 0.0f) {
        // T: Current world transformation of this joint
        glm::mat4 T = character.joint_nodes[joint_idx]->GetTransform().GetLocalToWorldMatrix();

        // Bi_inv: Inverse bind pose transformation
        glm::mat4 Bi_inv = glm::inverse(character.bind_pose_matrices[joint_idx]);

        // Ti = T * Bi^(-1) : The deformation transformation
        glm::mat4 Ti = T * Bi_inv;

        // p' = w * Ti * p
        deformed_pos += weight * (Ti * glm::vec4(original_pos, 1.0f));
        total_weight += weight;
      }
    }

    // 2. Handle Root Joint (Joint 0) with implicit weight
    float root_weight = 1.0f - total_weight;
    if (root_weight > 0.0f) {
      // T: Current world transformation of the root joint
      glm::mat4 T_root = character.joint_nodes[0]->GetTransform().GetLocalToWorldMatrix();

      // Bi_inv: Inverse bind pose transformation of the root
      glm::mat4 Bi_inv_root = glm::inverse(character.bind_pose_matrices[0]);

      // Ti = T * Bi^(-1) : The deformation transformation for the root
      glm::mat4 Ti_root = T_root * Bi_inv_root;

      // p' += w_root * Ti_root * p
      deformed_pos += root_weight * (Ti_root * glm::vec4(original_pos, 1.0f));
    }

    deformed_vertices[vertex_idx] = glm::vec3(deformed_pos);
  }

  // Compute normals
  std::vector<glm::vec3> deformed_normals(deformed_vertices.size());
  ComputeNormals(deformed_vertices, deformed_normals);

  // Update or create the deformed mesh object
  if (!character.deformed_mesh) {
    character.deformed_mesh = std::make_shared<VertexObject>();
    // Copy indices once
    if (bind_pose_mesh_->GetIndices().size() > 0) {
      character.deformed_mesh->UpdateIndices(GLOO::make_unique<IndexArray>(bind_pose_mesh_->GetIndices()));
    }
  }

  // Update positions and normals
  character.deformed_mesh->UpdatePositions(GLOO::make_unique<PositionArray>(deformed_vertices));
  character.deformed_mesh->UpdateNormals(GLOO::make_unique<NormalArray>(deformed_normals));

  // Set up the RenderingComponent for the SSD node if it hasn't been added yet
  if (character.ssd_node && !character.ssd_node->GetComponentPtr<RenderingComponent>()) {
    character.ssd_node->CreateComponent<RenderingComponent>(character.deformed_mesh);
  }
}

void SkeletonNode::LoadAttachmentWeights(const std::string& path) {
  std::ifstream file(path);
  std::string line;

  // Clear any existing weights
  attachment_weights_.clear();

  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::vector<float> row;
    float weight;

    // Parse each weight in the line
    while (iss >> weight) {
      row.push_back(weight);
    }

    // Only add non-empty rows
    if (!row.empty()) {
      attachment_weights_.push_back(row);
    }
  }

  file.close();
}

void SkeletonNode::InitializeCharacter(const std::string& prefix, const glm::vec3& world_pos) {
  // Create new character instance
  CharacterInstance character;

  // Create root node for world positioning
  auto root = make_unique<SceneNode>();
  root->GetTransform().SetPosition(world_pos);

  // Store raw pointer before moving
  SceneNode* root_ptr = root.get();
  AddChild(std::move(root));
  character.root_node.reset(root_ptr);

  // Load skeleton hierarchy
  std::string prefix_full = GetAssetDir() + prefix;
  LoadSkeletonFile(prefix_full + ".skel", character);

  // Compute bind pose matrices
  ComputeBindPoseMatrices(character);

  // Decorate tree with rendering nodes
  DecorateTree(character);

  // Initialize current pose with default values
  character.current_pose.resize(character.joint_nodes.size());
  for (auto& angle : character.current_pose) {
    angle.rx = angle.ry = angle.rz = 0.0f;
  }

  // Add to characters list
  characters_.push_back(std::move(character));
}

void SkeletonNode::LoadAllFiles(const std::string& prefix) {
  std::string prefix_full = GetAssetDir() + prefix;

  // Load shared resources only
  LoadMeshFile(prefix + ".obj");
  LoadAttachmentWeights(prefix_full + ".attach");

  // Setup shared rendering resources
  shader_ = std::make_shared<PhongShader>();
  float joint_radius = 0.03f;
  sphere_mesh_ = PrimitiveFactory::CreateSphere(joint_radius, 15, 15);
  float bone_radius = 0.015f;
  cylinder_mesh_ = PrimitiveFactory::CreateCylinder(bone_radius, 1.0f, 10);
}

void SkeletonNode::UpdateActiveCharacterPose(int new_idx) {
  if (new_idx < 0 || new_idx >= static_cast<int>(characters_.size())) {
    return;
  }

  // Save current UI values to the old character's pose
  if (active_char_idx_ >= 0 && active_char_idx_ < static_cast<int>(characters_.size())) {
    CharacterInstance& old_char = characters_[active_char_idx_];
    for (size_t i = 0; i < linked_angles_.size() && i < old_char.current_pose.size(); ++i) {
      old_char.current_pose[i] = *linked_angles_[i];
    }
  }

  // Switch active character
  active_char_idx_ = new_idx;

  // Load new character's pose into UI
  CharacterInstance& new_char = characters_[active_char_idx_];
  for (size_t i = 0; i < linked_angles_.size() && i < new_char.current_pose.size(); ++i) {
    *linked_angles_[i] = new_char.current_pose[i];
  }

  // Force scene update
  OnJointChanged();
}
}  // namespace GLOO
