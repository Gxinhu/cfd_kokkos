#pragma once
#include "Kokkos_Macros.hpp"
#include "util/cfd_shared.hpp"
#include "mesh/mesh_shared.hpp"
#include "util/parameters.hpp"
#include <cmath>
#include <memory>

namespace cfd_kokkos::mesh {

class NACA4DigitBase {
 public:
  explicit NACA4DigitBase(const MeshParams::MeshParamsPtr& mesh_params);
  KOKKOS_INLINE_FUNCTION
  precision calculate(precision x) const;
  void common_init();
  virtual void init() = 0;
  using MeshPtr       = std::shared_ptr<NACA4DigitBase>;
  size_t n_numbers_{};
  size_t m_numbers_{};
  MeshMatrix::HostMirror h_points_{};
  precision radius_{};
  precision chord_{};
  MeshMatrix points1_{};
  MeshMatrix points2_{};

 private:
  constexpr static precision kA0 = 0.2969;
  constexpr static precision kA1 = -0.1221;
  constexpr static precision kA2 = -0.3576;
  constexpr static precision kA3 = 0.2843;
  constexpr static precision kA4 = -0.1015;
  precision thickness_{};
  precision max_camber_{};
  precision max_location_{};
};

class NACA4DigitSys : public NACA4DigitBase {
 public:
  explicit NACA4DigitSys(const MeshParams::MeshParamsPtr& mesh_params);
  void init() override;
};

class NACA4DigitNonSys : public NACA4DigitBase {
 public:
  explicit NACA4DigitNonSys(const MeshParams::MeshParamsPtr& mesh_params);
  void init() override;
};

class MeshAirfoilFactory {
 public:
  MeshAirfoilFactory();
  MeshAirfoilFactory(const MeshAirfoilFactory&) = delete;
  static MeshAirfoilFactory& Instance();
  static std::shared_ptr<NACA4DigitBase> create(
      const MeshParams::MeshParamsPtr& mesh_params);
};

}  // namespace cfd_kokkos::mesh
