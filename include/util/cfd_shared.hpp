#pragma once
#include <json.hpp>

#include <Kokkos_Core.hpp>
#include "Kokkos_Macros.hpp"

namespace cfd_kokkos {
using Device = Kokkos::DefaultExecutionSpace;
using json   = nlohmann::json;

#ifdef USE_DOUBLE
using precision = double;
#else
using precision = float;
#endif

KOKKOS_INLINE_FUNCTION
void index2coord(int index, int &i, int &j, int Nx, int Ny) {
#ifdef KOKKOS_ENABLE_CUDA
  j = index / Nx;
  i = index - j * Nx;
#else
  i = index / Ny;
  j = index - i * Ny;
#endif
}

KOKKOS_INLINE_FUNCTION
int coord2index(int i, int j, int Nx, int Ny) {
#ifdef KOKKOS_ENABLE_CUDA
  return i + Nx * j;  // left layout
#else
  return j + Ny * i;  // right layout
#endif
}

/* 3D */

KOKKOS_INLINE_FUNCTION
void index2coord(int index, int &i, int &j, int &k, int Nx, int Ny, int Nz) {
#ifdef KOKKOS_ENABLE_CUDA
  int nx_ny = Nx * Ny;
  k         = index / nx_ny;
  j         = (index - k * nx_ny) / Nx;
  i         = index - j * Nx - k * nx_ny;
#else
  int NyNz = Ny * Nz;
  i        = index / NyNz;
  j        = (index - i * NyNz) / Nz;
  k        = index - j * Nz - i * NyNz;
#endif
}

KOKKOS_INLINE_FUNCTION
int coord2index(int i, int j, int k, int Nx, int Ny, int Nz) {
#ifdef KOKKOS_ENABLE_CUDA
  return i + Nx * j + Nx * Ny * k;  // left layout
#else
  return k + Nz * j + Nz * Ny * i;  // right layout
#endif
}

}  // namespace cfd_kokkos
