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
  void LoadAllFiles(const std::string& prefix);
  void LoadSkeletonFile(const std::string& path);
  void LoadMeshFile(const std::string& filename);
  void LoadAttachmentWeights(const std::string& path);
  void ComputeBindPoseMatrices();
  void DeformMesh();
  void ComputeNormals(const std::vector<glm::vec3>& vertices, std::vector<glm::vec3>& normals);

  void ToggleDrawMode();
  void DecorateTree();

  DrawMode draw_mode_;
  // Euler angles of the UI sliders.
  std::vector<EulerAngle*> linked_angles_;
  // Joint nodes for skeleton hierarchy.
  std::vector<SceneNode*> joint_nodes_;

  // Shared meshes and shader for stick figure rendering
  std::shared_ptr<VertexObject> sphere_mesh_;
  std::shared_ptr<VertexObject> cylinder_mesh_;
  std::shared_ptr<PhongShader> shader_;
  std::vector<SceneNode*> sphere_nodes_;
  std::vector<SceneNode*> cylinder_nodes_;

  // Bind pose mesh for SSD mode
  std::shared_ptr<VertexObject> bind_pose_mesh_;

  // Attachment weights matrix (n vertices × m-1 joints, excluding root joint)
  // attachment_weights_[i][j] = weight of joint (j+1) for vertex i
  // Root joint (joint 0) has implicit zero weights
  std::vector<std::vector<float>> attachment_weights_;

  // SSD deformation data
  std::vector<glm::mat4> bind_pose_matrices_;           // Bi matrices for each joint
  std::shared_ptr<VertexObject> deformed_mesh_;         // Deformed mesh for SSD mode
  SceneNode* ssd_node_;                                  // Node for SSD rendering

};
}  // namespace GLOO

#endif
