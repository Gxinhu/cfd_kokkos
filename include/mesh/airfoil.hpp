#pragma once
#include "Kokkos_Macros.hpp"
#include "util/cfd_shared.hpp"
#include "mesh/mesh_shared.hpp"
#include "util/parameters.hpp"
#include <cmath>

namespace cfd_kokkos::mesh {

class NACA4DigitBase {
 public:
  explicit NACA4DigitBase(const MeshParams::MeshPtr& mesh_params);
  KOKKOS_INLINE_FUNCTION
  precision calculate(precision x) const;
  void common_init();
  virtual void init() = 0;
  MeshMatrix points_{};

 protected:
  size_t numbers_{};
  MeshMatrix::HostMirror h_points_{};
  precision radius_{};
  precision chord_{};

 private:
  constexpr static precision kA0 = 0.2969;
  constexpr static precision kA1 = -0.1260;
  constexpr static precision kA2 = -0.3516;
  constexpr static precision kA3 = 0.2843;
  constexpr static precision kA4 = -0.1015;
  precision thickness_{};
  precision max_camber_{};
  precision max_location_{};
};

class NACA4DigitSys : public NACA4DigitBase {
 public:
  explicit NACA4DigitSys(const MeshParams::MeshPtr& mesh_params);
  void init() override;
};

class NACA4DigitNonSys : public NACA4DigitBase {
 public:
  explicit NACA4DigitNonSys(const MeshParams::MeshPtr& mesh_params);
  void init() override;
};

class MeshAirfoilFactory {
 public:
  MeshAirfoilFactory();
  MeshAirfoilFactory(const MeshAirfoilFactory&) = delete;
  static MeshAirfoilFactory& Instance();
  static std::shared_ptr<NACA4DigitBase> create(
      const MeshParams::MeshPtr& mesh_params);
};

}  // namespace cfd_kokkos::mesh
