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
  static std::shared_ptr<MeshSolverBase> create(
      const MeshParams::MeshParamsPtr &mesh_ptr);

 private:
  NodePtr new_point(FacePtr face, precision sp);
  NodePtr generate_tri(FacePtr face);

  /**
   * Finds candidate nodes for the front of a given node.
   *
   * @param new_node The node to find candidate nodes for.
   * @param sp The search precision.
   * @param face The face to search from.
   * @return A set of candidate nodes for the front of the given node.
   */
  std::set<NodePtr> find_candidate_node_face(const NodePtr &new_node,
                                             precision radius_pow,
                                             const FacePtr &face);
  /**
   * Finds candidate nodes for the front of the advancing front algorithm.
   *
   * @param new_node The new node being added to the mesh.
   * @param radius_pow The radius of the search area, raised to the power of 2.
   * @param face The face being processed.
   * @return A set of candidate nodes for the front of the advancing front
   * algorithm.
   */
  std::set<NodePtr> find_candidate_node_front(const NodePtr &new_node,
                                              precision radius_pow,
                                              const FacePtr &face);
  /**
   * Finds the candidate front of faces adjacent to the given set of candidate
   * nodes.
   *
   * @param candidates The set of candidate nodes to search for adjacent faces.
   * @return The set of candidate faces adjacent to the given set of candidate
   * nodes.
   */
  std::set<FacePtr> find_candidate_front(const std::set<NodePtr> &candidates);

  /**
   * @brief Finds the candidate faces for a given set of candidate nodes.
   *
   * @param candidates The set of candidate nodes.
   * @return The set of candidate faces.
   */
  std::set<FacePtr> find_candicate_face(const std::set<NodePtr> &candidates);

  std::set<CellPtr, CellComparator> quality_check_tri(
      const FacePtr &face, const std::set<NodePtr> &candidate_node);
  void update_info_tri(const FacePtr &face, const NodePtr &selected_node);

  void update_tri_cells(const FacePtr &face, const NodePtr &new_node);
  void remove_inactive_front();
  std::unique_ptr<Graph> graph_;
  precision search_range_ = 3.0;
  precision coeff_        = 0.8;
};

}  // namespace cfd_kokkos::mesh