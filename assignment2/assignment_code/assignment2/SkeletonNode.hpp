#ifndef SKELETON_NODE_H_
#define SKELETON_NODE_H_

#include "gloo/SceneNode.hpp"
#include "gloo/VertexObject.hpp"
#include "gloo/shaders/ShaderProgram.hpp"
#include "gloo/shaders/PhongShader.hpp"

#include <string>
#include <vector>
#include <memory>
#include <glm/glm.hpp>

namespace GLOO {
class SkeletonNode : public SceneNode {
 public:
  enum class DrawMode { Skeleton, SSD };
  struct EulerAngle {
    float rx, ry, rz;
  };

  SkeletonNode(const std::string& filename);
  void LinkRotationControl(const std::vector<EulerAngle*>& angles);
  void Update(double delta_time) override;
  void OnJointChanged();

 private:
  // Character instance - holds all per-character state
  struct CharacterInstance {
    std::unique_ptr<SceneNode> root_node;              // Root node for world positioning
    std::vector<SceneNode*> joint_nodes;               // Joint hierarchy
    std::vector<SceneNode*> sphere_nodes;              // Sphere rendering nodes
    std::vector<SceneNode*> cylinder_nodes;            // Cylinder rendering nodes
    SceneNode* ssd_node;                               // SSD mesh node
    std::shared_ptr<VertexObject> deformed_mesh;       // Deformed mesh
    std::vector<glm::mat4> bind_pose_matrices;         // Bi matrices
    std::vector<EulerAngle> current_pose;              // Current pose state
  };

  void LoadAllFiles(const std::string& prefix);
  void LoadSkeletonFile(const std::string& path, CharacterInstance& character);
  void LoadMeshFile(const std::string& filename);
  void LoadAttachmentWeights(const std::string& path);
  void ComputeBindPoseMatrices(CharacterInstance& character);
  void DeformMesh(CharacterInstance& character);
  void ComputeNormals(const std::vector<glm::vec3>& vertices, std::vector<glm::vec3>& normals);

  void ToggleDrawMode();
  void DecorateTree(CharacterInstance& character);
  void InitializeCharacter(const std::string& prefix, const glm::vec3& world_pos);
  void UpdateActiveCharacterPose(int new_idx);

  int active_char_idx_;
  DrawMode draw_mode_;
  // Euler angles of the UI sliders.
  std::vector<EulerAngle*> linked_angles_;

  // Multiple character instances
  std::vector<CharacterInstance> characters_;

  // Shared resources across all characters
  std::shared_ptr<VertexObject> sphere_mesh_;
  std::shared_ptr<VertexObject> cylinder_mesh_;
  std::shared_ptr<PhongShader> shader_;
  std::shared_ptr<VertexObject> bind_pose_mesh_;
  std::vector<std::vector<float>> attachment_weights_;

};
}  // namespace GLOO

#endif
