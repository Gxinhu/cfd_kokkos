#pragma once

#include "util/cfd_shared.hpp"
#include <Kokkos_Core.hpp>
namespace cfd_kokkos::mesh {
using MeshPoints = Kokkos::View<precision **, Device>;
using MeshMatrix = Kokkos::View<precision ***, Device>;
}  // namespace cfd_kokkos::mesh
