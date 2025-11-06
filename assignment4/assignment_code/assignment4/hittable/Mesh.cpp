#include "Mesh.hpp"

#include <functional>
#include <stdexcept>
#include <iostream>

#include "gloo/utils.hpp"

namespace GLOO {
Mesh::Mesh(std::unique_ptr<PositionArray> positions,
           std::unique_ptr<NormalArray> normals,
           std::unique_ptr<IndexArray> indices) {
  size_t num_vertices = indices->size();
  if (num_vertices % 3 != 0 || normals->size() != positions->size())
    throw std::runtime_error("Bad mesh data in Mesh constuctor!");

  for (size_t i = 0; i < num_vertices; i += 3) {
    triangles_.emplace_back(
        positions->at(indices->at(i)), positions->at(indices->at(i + 1)),
        positions->at(indices->at(i + 2)), normals->at(indices->at(i)),
        normals->at(indices->at(i + 1)), normals->at(indices->at(i + 2)));
  }
  // Let mesh data destruct.

  // Build Octree.
  octree_ = make_unique<Octree>();
  octree_->Build(*this);
}

Mesh::Mesh(std::unique_ptr<PositionArray> positions,
           std::unique_ptr<NormalArray> normals,
           std::unique_ptr<IndexArray> indices,
           std::unique_ptr<TexCoordArray> tex_coords) {
  size_t num_vertices = indices->size();
  if (num_vertices % 3 != 0 || normals->size() != positions->size())
    throw std::runtime_error("Bad mesh data in Mesh constuctor!");

  bool has_tex_coords = (tex_coords != nullptr && tex_coords->size() == positions->size());

  for (size_t i = 0; i < num_vertices; i += 3) {
    size_t idx0 = indices->at(i);
    size_t idx1 = indices->at(i + 1);
    size_t idx2 = indices->at(i + 2);

    if (has_tex_coords) {
      // Create triangle with texture coordinates
      triangles_.emplace_back(
          positions->at(idx0), positions->at(idx1), positions->at(idx2),
          normals->at(idx0), normals->at(idx1), normals->at(idx2),
          tex_coords->at(idx0), tex_coords->at(idx1), tex_coords->at(idx2));
    } else {
      // Create triangle without texture coordinates (uses default UVs)
      triangles_.emplace_back(
          positions->at(idx0), positions->at(idx1), positions->at(idx2),
          normals->at(idx0), normals->at(idx1), normals->at(idx2));
    }
  }
  // Let mesh data destruct.

  // Build Octree.
  octree_ = make_unique<Octree>();
  octree_->Build(*this);
}

bool Mesh::Intersect(const Ray& ray, float t_min, HitRecord& record) const {
  return octree_->Intersect(ray, t_min, record);
}
}  // namespace GLOO
