#pragma once

#include <Kokkos_Core.hpp>
#include <utility>

#include "util/cfd_shared.hpp"
namespace cfd_kokkos::mesh {
using MeshPoints = Kokkos::View<precision **, Device>;
using MeshMatrix = Kokkos::View<precision ***, Device>;
enum BoundaryType { kRealBoundary, kNewBoundary };
}  // namespace cfd_kokkos::mesh
