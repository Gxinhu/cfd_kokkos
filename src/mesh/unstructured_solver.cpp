#include <algorithm>
#include <cstddef>
#include <iterator>
#include <memory>
#include <set>

#include "fmt/core.h"
#include "fmt/format.h"
#include "io/save_vtk.hpp"
#include "mesh/common_function.hpp"
#include "mesh/functor.hpp"
#include "mesh/graph.hpp"
#include "mesh/mesh_shared.hpp"
#include "mesh/solver.hpp"
#include "mesh/unsturctured_solver.hpp"
#include "util/cfd_shared.hpp"
#include "util/parameters.hpp"

namespace cfd_kokkos::mesh {
std::shared_ptr<MeshSolverBase> MeshSolverOUnstructured::create(
    const MeshParams::MeshParamsPtr& mesh_ptr) {
  auto solver = std::make_shared<MeshSolverOUnstructured>(mesh_ptr);
  return solver;
}

MeshSolverOUnstructured::MeshSolverOUnstructured(
    const MeshParams::MeshParamsPtr& mesh_ptr)
    : MeshSolverO(mesh_ptr) {
  graph_ = std::make_unique<Graph>();
}
void MeshSolverOUnstructured::init() {
  fmt::print("###Init Mesh###");

  precision d_theta = (M_PI * 2) / mesh_->n_numbers_;
  for (int i = 0; i < mesh_->n_numbers_; ++i) {
    auto node = std::make_shared<Node>(i, mesh_->h_points_(0, i, 0),
                                       mesh_->h_points_(0, i, 1));
    graph_->add_node(node);
    graph_->x_min_ = std::min(graph_->x_min_, node->x_);
    graph_->x_max_ = std::max(graph_->x_max_, node->x_);
    graph_->y_min_ = std::min(graph_->y_min_, node->y_);
    graph_->y_max_ = std::max(graph_->y_max_, node->y_);
  }
  for (int i = 1; i < graph_->node_num_; ++i) {
    auto face = std::make_shared<Face>(graph_->nodes_[i - 1], graph_->nodes_[i],
                                       BoundaryType::kRealBoundary, -1, 0);
    graph_->add_front(face);
  }
  auto last_node        = graph_->nodes_[graph_->node_num_ - 1];
  auto surface_node_num = graph_->node_num_;
  graph_->add_front(std::make_shared<Face>(last_node, graph_->nodes_[0],
                                           BoundaryType::kRealBoundary, -1, 0));

  for (int i = 0; i < mesh_->n_numbers_; ++i) {
    auto theta = i * d_theta;
    auto x     = mesh_->chord_ / 2.0 + mesh_->radius_ * std::cos(theta);
    auto y     = -mesh_->radius_ * std::sin(theta);
    auto node  = std::make_shared<Node>(i + surface_node_num, x, y);
    graph_->add_node(node);

    graph_->x_min_ = std::min(graph_->x_min_, node->x_);
    graph_->x_max_ = std::max(graph_->x_max_, node->x_);
    graph_->y_min_ = std::min(graph_->y_min_, node->y_);
    graph_->y_max_ = std::max(graph_->y_max_, node->y_);
  }
  last_node = graph_->nodes_[surface_node_num];

  for (int i = surface_node_num + 1; i < graph_->node_num_; ++i) {
    auto node = graph_->nodes_[i];
    auto face = std::make_shared<Face>(last_node, node,
                                       BoundaryType::kRealBoundary, -1, 0);
    graph_->add_front(face);
    last_node = node;
  }
  graph_->add_front(std::make_shared<Face>(last_node,
                                           graph_->nodes_[surface_node_num],
                                           BoundaryType::kRealBoundary, -1, 0));
}

void MeshSolverOUnstructured::solve() {
  while (graph_->front_num_ > 0) {
    auto selected_node = generate_tri(*graph_->front_list_.begin());
    if (selected_node == nullptr) {
      search_range_ = search_range_ * 1.2;
    } else {
      search_range_ = 3.0;
    }
  }
  // TODO(xhu): 在这里实现阵面推进算法
}

NodePtr MeshSolverOUnstructured::new_point(FacePtr face) {
  auto normal_vector_x =
      -(face->node2_->y_ - face->node1_->y_) / (face->distance_ + 1e-10);

  auto normal_vector_y =
      -(face->node2_->x_ - face->node1_->x_) / (face->distance_ + 1e-10);

  auto x        = (face->node1_->x_ + face->node2_->x_) / 2.0 + normal_vector_x;
  auto y        = (face->node1_->y_ + face->node2_->y_) / 2.0 + normal_vector_y;
  auto new_node = std::make_shared<Node>(-1, x, y);
  return new_node;
}

NodePtr MeshSolverOUnstructured::generate_tri(FacePtr face) {
  auto sp = face->distance_ * 0.866;
  if (face->boundary_type_ == BoundaryType::kRealBoundary) {
    sp = face->distance_ * sqrt(3.0) / 2.0;
  }
  auto new_node   = new_point(face);
  auto radius_pow = pow(search_range_ * sp, 2);
  auto candidates_node_front =
      find_candicate_node_front(new_node, radius_pow, face);
  auto candidate_front = find_candicate_front(candidates_node_front);

  auto candidates_node_face =
      find_candicate_node_face(new_node, radius_pow, face);
  std::set<NodePtr> candidate_node;
  std::set_union(candidates_node_front.begin(), candidates_node_front.end(),
                 candidates_node_face.begin(), candidates_node_face.end(),
                 std::inserter(candidate_node, candidate_node.begin()));
  auto candidates_face  = find_candicate_face(candidate_node);
  auto candidate_cell   = quality_check_tri(face, candidates_node_front);
  CellPtr selected_cell = nullptr;
  for (const auto& cell : candidate_cell) {
    bool flag = false;
    for (const auto& front : candidate_front) {
      if (is_cross(cell->node1_, cell->node3_, front->node1_, front->node2_)) {
        flag = true;
        break;
      }
      if (is_cross(cell->node1_, cell->node2_, front->node1_, front->node2_)) {
        flag = true;
        break;
      }
    }
    if (flag) {
      continue;
    }

    for (const auto& front : candidates_face) {
      if (is_cross(cell->node1_, cell->node3_, front->node1_, front->node2_)) {
        flag = true;
        break;
      }
      if (is_cross(cell->node1_, cell->node2_, front->node1_, front->node2_)) {
        flag = true;
        break;
      }
    }
    if (flag) {
      continue;
    }
    if (cell->is_right_cell_) {
      selected_cell = cell;
      // FIXME 目前左右 cell 还有问题，需要解决
    }
  }
  if (selected_cell != nullptr) {
    if (selected_cell->node3_->id_ == -1) {
      graph_->add_node(selected_cell->node3_);
    }
    return selected_cell->node3_;
  }
  return nullptr;
}
std::set<NodePtr> MeshSolverOUnstructured::find_candicate_node_front(
    const NodePtr& new_node, precision radius_pow, const FacePtr& face) {
  std::set<NodePtr> candidates;
  candidates.insert(new_node);

  for (const auto& front : graph_->front_list_) {
    auto distance1 = pow(front->node1_->x_ - new_node->x_, 2) +
                     pow(front->node1_->y_ - new_node->y_, 2);
    auto distance2 = pow(front->node2_->x_ - new_node->x_, 2) +
                     pow(front->node2_->y_ - new_node->y_, 2);
    if (distance1 < radius_pow && front->node2_ != face->node1_ &&
        front->node1_ != face->node1_) {
      candidates.insert(front->node1_);
    }
    if (distance2 < radius_pow && front->node2_ != face->node2_ &&
        front->node1_ != face->node2_) {
      candidates.insert(front->node2_);
    }
    if ((front->node2_ == face->node2_ && front->node1_ != face->node1_) ||
        (front->node2_ == face->node1_ && front->node1_ != face->node2_)) {
      candidates.insert(front->node1_);
    }
    if ((front->node1_ == face->node2_ && front->node2_ != face->node1_) ||
        (front->node1_ == face->node1_ && front->node2_ != face->node2_)) {
      candidates.insert(front->node2_);
    }
  }
  return candidates;
}

std::set<NodePtr> MeshSolverOUnstructured::find_candicate_node_face(
    const NodePtr& new_node, precision radius_pow, const FacePtr& face) {
  std::set<NodePtr> candidates;
  candidates.insert(new_node);

  for (const auto& front : graph_->faces_) {
    auto distance1 = pow(front->node1_->x_ - new_node->x_, 2) +
                     pow(front->node1_->y_ - new_node->y_, 2);
    auto distance2 = pow(front->node2_->x_ - new_node->x_, 2) +
                     pow(front->node2_->y_ - new_node->y_, 2);
    if (distance1 < radius_pow && front->node2_ != face->node1_ &&
        front->node1_ != face->node1_) {
      candidates.insert(front->node1_);
    }
    if (distance2 < radius_pow && front->node2_ != face->node2_ &&
        front->node1_ != face->node2_) {
      candidates.insert(front->node2_);
    }
    if ((front->node2_ == face->node2_ && front->node1_ != face->node1_) ||
        (front->node2_ == face->node1_ && front->node1_ != face->node2_)) {
      candidates.insert(front->node1_);
    }
    if ((front->node1_ == face->node2_ && front->node2_ != face->node1_) ||
        (front->node1_ == face->node1_ && front->node2_ != face->node2_)) {
      candidates.insert(front->node2_);
    }
  }
  return candidates;
}

std::set<FacePtr> MeshSolverOUnstructured::find_candicate_front(
    const std::set<NodePtr>& candidates) {
  std::set<FacePtr> candidate_front;
  for (auto& front : graph_->front_list_) {
    if (candidates.find(front->node1_) != candidates.end() ||
        candidates.find(front->node2_) != candidates.end()) {
      candidate_front.insert(front);
    }
  }
  return candidate_front;
}

std::set<FacePtr> MeshSolverOUnstructured::find_candicate_face(
    const std::set<NodePtr>& candidates) {
  std::set<FacePtr> candidate_front;
  for (auto& front : graph_->faces_) {
    if (candidates.find(front->node1_) != candidates.end() ||
        candidates.find(front->node2_) != candidates.end()) {
      candidate_front.insert(front);
    }
  }
  return candidate_front;
}
std::set<CellPtr, CellComparator> MeshSolverOUnstructured::quality_check_tri(
    const FacePtr& face, const std::set<NodePtr>& candidate_node) {
  std::set<CellPtr, CellComparator> candidate_cell;
  for (auto& node : candidate_node) {
    auto cell = std::make_shared<Cell>(-1, face->node1_, face->node2_, node);
    if (node->id_ == -1) {
      cell->quality_ = coeff_;
    }
    candidate_cell.insert(cell);
  }
  return candidate_cell;
}

}  // namespace cfd_kokkos::mesh