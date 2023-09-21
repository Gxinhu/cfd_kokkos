#include "mesh/graph.hpp"

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

}  // namespace cfd_kokkos::mesh