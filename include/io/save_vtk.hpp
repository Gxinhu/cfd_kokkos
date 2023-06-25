#pragma once

#include "mesh/airfoil.hpp"
#include "util/cfd_shared.hpp"
#include "mesh/mesh_shared.hpp"
#include "util/parameters.hpp"
namespace cfd_kokkos::io {
void save_vtk_2d(const mesh::MeshMatrix::HostMirror &h_points,
                 const std::string &save_name);
}