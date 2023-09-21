
#pragma once

#include "mesh/graph.hpp"
namespace cfd_kokkos::mesh {
bool is_cross(const NodePtr &a, const NodePtr &b, const NodePtr &c,
              const NodePtr &d);

}  // namespace cfd_kokkos::mesh