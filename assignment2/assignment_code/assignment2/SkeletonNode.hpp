#ifndef SKELETON_NODE_H_
#define SKELETON_NODE_H_

#include "gloo/SceneNode.hpp"
#include "gloo/VertexObject.hpp"
#include "gloo/shaders/ShaderProgram.hpp"

#include <string>
#include <vector>

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

  void ToggleDrawMode();
  void DecorateTree();

  DrawMode draw_mode_;
  // Euler angles of the UI sliders.
  std::vector<EulerAngle*> linked_angles_;
};
}  // namespace GLOO

#endif
