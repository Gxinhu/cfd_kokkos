#pragma once

#include <memory>

#include "mesh/graph.hpp"
#include "mesh/solver.hpp"
#include "util/cfd_shared.hpp"

namespace cfd_kokkos::mesh {

class MeshSolverOUnstructured : public MeshSolverO {
 public:
  explicit MeshSolverOUnstructured(const MeshParams::MeshParamsPtr &mesh_ptr);
  void solve() override;
  void init() override;
  NodePtr new_point(FacePtr face);
  NodePtr generate_tri(FacePtr face);

  /**
   * Finds candidate nodes for the front of a given node.
   *
   * @param new_node The node to find candidate nodes for.
   * @param sp The search precision.
   * @param face The face to search from.
   * @return A set of candidate nodes for the front of the given node.
   */
  std::set<NodePtr> find_candicate_node_face(const NodePtr &new_node,
                                             precision radius_pow,
                                             const FacePtr &face);
  std::set<NodePtr> find_candicate_node_front(const NodePtr &new_node,
                                              precision radius_pow,
                                              const FacePtr &face);
  std::set<FacePtr> find_candicate_front(const std::set<NodePtr> &candidates);

  std::set<FacePtr> find_candicate_face(const std::set<NodePtr> &candidates);

  std::set<CellPtr, CellComparator> quality_check_tri(
      const FacePtr &face, const std::set<NodePtr> &candidate_node);
  static std::shared_ptr<MeshSolverBase> create(
      const MeshParams::MeshParamsPtr &mesh_ptr);

 private:
  std::unique_ptr<Graph> graph_;
  precision search_range_ = 3.0;
  precision coeff_        = 0.8;
};

}  // namespace cfd_kokkos::mesh