#include "mesh/common_function.hpp"

namespace cfd_kokkos::mesh {

bool is_cross(const NodePtr &a, const NodePtr &b, const NodePtr &c,
              const NodePtr &d) {
  if ((std::min(a->x_, b->x_) <= std::max(c->x_, d->x_)) and
      (std::min(c->x_, d->x_) <= std::max(a->x_, b->x_)) and
      (std::min(a->y_, b->y_) <= std::max(c->y_, d->y_)) and
      (std::min(c->y_, d->y_) <= std::max(a->y_, b->y_))) {
    auto u =
        (c->x_ - a->x_) * (b->y_ - a->y_) - (b->x_ - a->x_) * (c->y_ - a->y_);
    auto v =
        (d->x_ - a->x_) * (b->y_ - a->y_) - (b->x_ - a->x_) * (d->y_ - a->y_);
    auto w =
        (a->x_ - c->x_) * (d->y_ - c->y_) - (d->x_ - c->x_) * (a->y_ - c->y_);
    auto z =
        (b->x_ - c->x_) * (d->y_ - c->y_) - (d->x_ - c->x_) * (b->y_ - c->y_);
    const auto eps = 1e-9;
    if ((u * v <= 0) && (w * z <= 0)) {
      return ((abs(a->x_ - c->x_) >= eps) || (abs(a->x_ - c->x_) >= eps)) &&
             ((abs(a->x_ - d->x_) >= eps) || (abs(a->x_ - d->x_) >= eps)) &&
             ((abs(b->x_ - c->x_) >= eps) || (abs(b->x_ - c->x_) >= eps)) &&
             ((abs(b->x_ - d->x_) >= eps) || (abs(b->x_ - d->x_) >= eps));
    }
  }
  return false;
}

}  // namespace cfd_kokkos::mesh