#pragma once

#include <limits>
#include <memory>
#include <set>
#include <unordered_set>
#include <utility>

#include "mesh/mesh_shared.hpp"
#include "mesh/solver.hpp"
#include "util/cfd_shared.hpp"

namespace cfd_kokkos::mesh {

class Node;
class Face;
class Cell;
using NodePtr = std::shared_ptr<Node>;
using FacePtr = std::shared_ptr<Face>;
using CellPtr = std::shared_ptr<Cell>;

class Node {
 public:
  explicit Node(int id, precision x, precision y) : id_(id), x_(x), y_(y){};
  int id_;
  precision x_;
  precision y_;
  std::set<NodePtr> neighbors_;
};

class Face {
 public:
  Face(NodePtr node1, NodePtr node2, BoundaryType boundary_type, int left_cell,
       int right_cell)
      : node1_(std::move(node1)),
        node2_(std::move(node2)),
        boundary_type_(boundary_type),
        left_cell_(left_cell),
        right_cell_(right_cell) {
    distance_ = std::sqrt(std::pow(node1_->x_ - node2_->x_, 2) +
                          std::pow(node1_->y_ - node2_->y_, 2));
    node1_->neighbors_.insert(node2_);
    node2_->neighbors_.insert(node1_);
  }
  precision distance_;
  int boundary_type_;
  int left_cell_;
  int right_cell_;
  NodePtr node1_;
  NodePtr node2_;
  NodePtr target_node_;
  bool operator<(const Face &other) const {
    return distance_ < other.distance_;
  }
};
class Cell {
 public:
  Cell(int id, NodePtr node1, NodePtr node2, NodePtr node3)
      : id_(id),
        node1_(std::move(node1)),
        node2_(std::move(node2)),
        node3_(std::move(node3)) {
    quality_check();
    is_leftcell();
  };
  void quality_check() {
    auto a   = sqrt(pow((node1_->x_ - node2_->x_), 2) +
                    pow((node1_->y_ - node2_->y_), 2));
    auto b   = sqrt(pow((node2_->x_ - node3_->x_), 2) +
                    pow((node2_->y_ - node3_->y_), 2));
    auto c   = sqrt(pow((node1_->x_ - node3_->x_), 2) +
                    pow((node1_->y_ - node3_->y_), 2));
    auto tmp = (pow(a, 2) + pow(b, 2) - pow(c, 2)) / (2.0 * a * b);
    if (abs(tmp - 1.0) < 1e-5) {
      tmp = 1;
    }
    if (abs(tmp + 1.0) < 1e-5) {
      tmp = -1;
    }
    auto theta = acos(tmp);
    auto area  = 0.5 * a * b * sin(theta) + 1e-40;
    quality_   = (4.0 * sqrt(3.0) * area / (pow(a, 2) + pow(b, 2) + pow(c, 2)));
  }
  void is_leftcell() {
    auto x1       = node2_->x_ - node1_->x_;
    auto y1       = node2_->y_ - node1_->y_;
    auto x2       = node3_->x_ - node1_->x_;
    auto y2       = node3_->y_ - node1_->y_;
    auto res      = x1 * y2 - x2 * y1;
    is_right_cell_ = (res > 0);
  }

  int id_;
  precision quality_;
  NodePtr node1_;
  NodePtr node2_;
  NodePtr node3_;
  bool is_right_cell_;

  bool operator>(const Cell &other) const { return quality_ > other.quality_; }
};

struct CellComparator {
  bool operator()(const CellPtr &a, const CellPtr &b) const { return a > b; }
};
struct FaceComparator {
  bool operator()(const FacePtr &a, const FacePtr &b) const { return a < b; }
};
class Graph {
 public:
  Graph() = default;
  std::vector<NodePtr> nodes_;
  std::set<FacePtr> faces_;
  std::set<CellPtr, CellComparator> cells_;
  std::set<FacePtr, FaceComparator> front_list_;
  int node_num_          = 0;
  int face_num_          = 0;
  int cell_num_          = 0;
  int front_num_         = 0;
  int front_new_         = 0;
  precision x_min_       = std::numeric_limits<precision>::max();
  precision y_min_       = std::numeric_limits<precision>::max();
  precision x_max_       = std::numeric_limits<precision>::min();
  precision y_max_       = std::numeric_limits<precision>::min();
  precision search_range = 3.0;

  void add_node(const NodePtr &node);
  void add_face(const FacePtr &face);
  void add_cell(const CellPtr &cell);
  void remove_cell(const CellPtr &cell);
  void add_front(const FacePtr &front);
  void remove_front(const FacePtr &front);
  void sort_front();

 private:
};
}  // namespace cfd_kokkos::mesh