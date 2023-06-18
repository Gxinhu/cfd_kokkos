#include "mesh/airfoil.hpp"
#include "mesh/mesh_shared.hpp"

namespace cfd_kokkos::mesh {

NACA4DigitBase::NACA4DigitBase(const MeshParams::MeshPtr& mesh_params)
    : chord_(mesh_params->chord_),
      thickness_(mesh_params->thickness_ / 100),
      max_camber_(mesh_params->max_camber_),
      max_location_(mesh_params->max_location_),
      numbers_(mesh_params->n_numbers_),
      radius_(mesh_params->radius_),
      points_(MeshMatrix("mesh_points", mesh_params->m_numbers_,
                         mesh_params->n_numbers_, 2)) {
  h_points_ = Kokkos::create_mirror(points_);
}

KOKKOS_INLINE_FUNCTION
precision NACA4DigitBase::calculate(precision x) const {
  return 5 * thickness_ * chord_ *
         (kA0 * std::sqrt(x) + kA1 * x + kA2 * std::pow(x, 2) +
          kA3 * std::pow(x, 3) + kA4 * std::pow(x, 4));
}

void NACA4DigitBase::common_init() {
  precision d_theta = (M_PI * 2) / (numbers_ - 1);
  for (int i = 0; i < numbers_ / 2; ++i) {
    auto theta         = i * d_theta;
    auto x             = (1.0 + std::cos(theta)) / 2.0 - chord_ / 2.0;
    auto y             = calculate(x);
    h_points_(0, i, 0) = x;
    h_points_(0, i, 1) = y;
  }
}

NACA4DigitSys::NACA4DigitSys(const MeshParams::MeshPtr& mesh_params)
    : NACA4DigitBase(mesh_params) {}

void NACA4DigitSys::init() {
  common_init();

  precision d_theta = (M_PI * 2) / (numbers_ - 1);
  for (int i = numbers_ / 2.0 + 1; i < numbers_; ++i) {
    auto theta         = i * d_theta;
    auto x             = (1.0 + std::cos(theta)) / 2.0 - chord_ / 2.0;
    auto y             = -calculate(x);
    h_points_(0, i, 0) = x;
    h_points_(0, i, 1) = y;
  }

  for (int i = 0; i < numbers_; ++i) {
    auto theta                    = i * d_theta;
    auto x                        = radius_ * std::cos(theta);
    auto y                        = -radius_ * std::sin(theta);
    h_points_(numbers_ - 1, i, 0) = x;
    h_points_(numbers_ - 1, i, 1) = y;
  }

  Kokkos::deep_copy(points_, h_points_);
}

NACA4DigitNonSys::NACA4DigitNonSys(const MeshParams::MeshPtr& mesh_params)
    : NACA4DigitBase(mesh_params) {}

void NACA4DigitNonSys::init() {
  common_init();
  // TODO(xhu): Implement Non Sys Naca airfoil code
  Kokkos::deep_copy(points_, h_points_);
}

MeshAirfoilFactory::MeshAirfoilFactory() = default;

MeshAirfoilFactory& MeshAirfoilFactory::Instance() {
  static MeshAirfoilFactory instance;
  return instance;
}

std::shared_ptr<NACA4DigitBase> MeshAirfoilFactory::create(
    const MeshParams::MeshPtr& mesh_params) {
  std::shared_ptr<NACA4DigitBase> mesh_ptr;
  if (mesh_params->is_sys_) {
    mesh_ptr = std::make_shared<NACA4DigitSys>(mesh_params);
  } else {
    mesh_ptr = std::make_shared<NACA4DigitNonSys>(mesh_params);
  }
  mesh_ptr->init();
  return mesh_ptr;
}

}  // namespace cfd_kokkos::mesh
