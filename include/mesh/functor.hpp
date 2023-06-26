
#include <utility>

#include "Kokkos_Macros.hpp"
#include "impl/Kokkos_Profiling_C_Interface.h"
#include "mesh/mesh_shared.hpp"
#include "util/cfd_shared.hpp"
#include "util/parameters.hpp"

namespace cfd_kokkos::mesh {

class BoundaryFitFunctor {
  BoundaryFitFunctor(MeshMatrix in_points, MeshMatrix out_points)
      : in_points_(std::move(in_points)), out_points_(std::move(out_points)){};

 public:
  static void apply(const MeshMatrix& in_points, const MeshMatrix& out_points,
                    precision& error) {
    BoundaryFitFunctor functor(in_points, out_points);
    Kokkos::parallel_reduce(in_points.extent(0) * in_points.extent(1), functor,
                            Kokkos::Max<precision>(error));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& index, precision& error) const {
    const auto nx = in_points_.extent(0);
    const auto ny = in_points_.extent(1);
    int i         = -1;
    int j         = -1;
    index2coord(index, i, j, nx, ny);
    if (i == 0 or i == nx - 1) {
      return;
    }

    precision x_right      = -1;
    precision x_left       = -1;
    precision x_down       = -1;
    precision x_up         = -1;
    precision x_up_right   = -1;
    precision x_up_left    = -1;
    precision x_down_right = -1;
    precision x_down_left  = -1;

    precision y_right      = -1;
    precision y_left       = -1;
    precision y_down       = -1;
    precision y_up         = -1;
    precision y_up_right   = -1;
    precision y_up_left    = -1;
    precision y_down_right = -1;
    precision y_down_left  = -1;
    if (j == 0) {
      x_right      = in_points_(i, j + 1, 0);
      x_up_right   = in_points_(i + 1, j + 1, 0);
      x_down_right = in_points_(i - 1, j + 1, 0);
      x_left       = in_points_(i, ny - 2, 0);
      x_up_left    = in_points_(i + 1, ny - 2, 0);
      x_down_left  = in_points_(i - 1, ny - 2, 0);
      y_right      = in_points_(i, j + 1, 1);
      y_up_right   = in_points_(i + 1, j + 1, 1);
      y_down_right = in_points_(i - 1, j + 1, 1);
      y_left       = in_points_(i, ny - 2, 1);
      y_up_left    = in_points_(i + 1, ny - 2, 1);
      y_down_left  = in_points_(i - 1, ny - 2, 1);
    } else if (j == ny - 1) {
      x_right      = in_points_(i, 1, 0);
      x_up_right   = in_points_(i + 1, 1, 0);
      x_down_right = in_points_(i - 1, 1, 0);
      x_left       = in_points_(i, j - 1, 0);
      x_up_left    = in_points_(i + 1, j - 1, 0);
      x_down_left  = in_points_(i - 1, j - 1, 0);
      y_right      = in_points_(i, 1, 1);
      y_up_right   = in_points_(i + 1, 1, 1);
      y_down_right = in_points_(i - 1, 1, 1);
      y_left       = in_points_(i, j - 1, 1);
      y_up_left    = in_points_(i + 1, j - 1, 1);
      y_down_left  = in_points_(i - 1, j - 1, 1);
    } else {
      x_right      = in_points_(i, j + 1, 0);
      x_left       = in_points_(i, j - 1, 0);
      x_up_right   = in_points_(i + 1, j + 1, 0);
      x_up_left    = in_points_(i + 1, j - 1, 0);
      x_down_right = in_points_(i - 1, j + 1, 0);
      x_down_left  = in_points_(i - 1, j - 1, 0);
      y_right      = in_points_(i, j + 1, 1);
      y_left       = in_points_(i, j - 1, 1);
      y_up_right   = in_points_(i + 1, j + 1, 1);
      y_up_left    = in_points_(i + 1, j - 1, 1);
      y_down_right = in_points_(i - 1, j + 1, 1);
      y_down_left  = in_points_(i - 1, j - 1, 1);
    }
    x_down     = in_points_(i - 1, j, 0);
    x_up       = in_points_(i + 1, j, 0);
    y_down     = in_points_(i - 1, j, 1);
    y_up       = in_points_(i + 1, j, 1);
    auto x_eta = x_up - x_down;
    auto x_xi  = x_right - x_left;
    auto y_eta = y_up - y_down;
    auto y_xi  = y_right - y_left;
    auto alpha = x_eta * x_eta / 4.0 + y_eta * y_eta / 4.0;
    auto beta  = x_eta * x_xi / 4.0 + y_eta * y_xi / 4.0;
    auto gamma = x_xi * x_xi / 4.0 + y_xi * y_xi / 4.0;

    auto xx_sum = x_right + x_left;
    auto xy_sum = x_up + x_down;
    auto yx_sum = y_right + y_left;
    auto yy_sum = y_up + y_down;

    auto x =
        0.5 *
        (alpha * xx_sum + gamma * xy_sum -
         0.5 * beta * (x_up_right + x_down_left - x_down_right - x_up_left)) /
        (alpha + gamma);

    auto y =
        0.5 *
        (alpha * yx_sum + gamma * yy_sum -
         0.5 * beta * (y_up_right + y_down_left - y_down_right - y_up_left)) /
        (alpha + gamma);
    auto x_error         = Kokkos::fabs(out_points_(i, j, 0) - x);
    auto y_error         = Kokkos::fabs(out_points_(i, j, 1) - y);
    error                = Kokkos::fmax(x_error, y_error);
    out_points_(i, j, 0) = x;
    out_points_(i, j, 1) = y;
  }

 private:
  MeshMatrix in_points_;
  MeshMatrix out_points_;
};

}  // namespace cfd_kokkos::mesh