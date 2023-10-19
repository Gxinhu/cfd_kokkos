#include "mesh/graph.hpp"

#include <utility>

namespace cfd_kokkos::mesh {

void Graph::add_node(const NodePtr &node) {
  node->id_ = node_num_;
  nodes_.push_back(node);
  ++node_num_;
}

void Graph::add_face(const FacePtr &face) {
  faces_.insert(face);
  ++face_num_;
}

void Graph::add_cell(const CellPtr &cell) {
  cells_.insert(cell);
  ++cell_num_;
}

void Graph::remove_cell(const CellPtr &cell) {
  cells_.erase(cell);
  --cell_num_;
}

void Graph::add_front(const FacePtr &front) {
  front_list_.insert(front);
  ++front_num_;
}

void Graph::remove_front(const FacePtr &front) {
  front_list_.erase(front);
  --front_num_;
}

std::pair<FacePtr, bool> Graph::front_exist(const NodePtr &node1,
                                            const NodePtr &node2) {
  FacePtr exist_face = nullptr;
  auto direction     = false;
  for (const auto &front : front_list_) {
    if (front->node1_ == node1 && front->node2_ == node2) {
      exist_face = front;
      direction  = true;
      break;
    }
    if (front->node2_ == node1 && front->node1_ == node2) {
      exist_face = front;
      direction  = false;
      break;
    }
  }
  return std::make_pair(exist_face, direction);
}
void Graph::remove_inactive_front() {
  auto iter=front_list_.begin();
  while (iter != front_list_.end()) {
    if ((*iter)->left_cell_ != -1 && (*iter)->right_cell_ != -1) {
      add_face(*iter);
      auto tmp = iter;
      ++iter;
      remove_front(*tmp);
    } else {
      ++iter;
    }
  }
};

}  // namespace cfd_kokkos::mesh