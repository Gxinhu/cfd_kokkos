#include "mesh/solver.hpp"
#include "io/save_vtk.hpp"

namespace cfd_kokkos::mesh {

MeshSolverBase::MeshSolverBase(const MeshParams::MeshParamsPtr& mesh_ptr) {
  mesh_ = cfd_kokkos::mesh::MeshAirfoilFactory::create(mesh_ptr);
}

void MeshSolverBase::save_mesh() const {
  io::save_vtk_2d(mesh_->h_points_, "mesh.vts");
}

MeshSolverO::MeshSolverO(const MeshParams::MeshParamsPtr& mesh_ptr)
    : MeshSolverBase(mesh_ptr) {}

void MeshSolverO::solve() { fmt::print("Hello World\n"); }

void MeshSolverO::init() {
  fmt::print("###Init Mesh###");

  precision d_theta = (M_PI * 2) / (mesh_->n_numbers_ - 1);
  for (int i = 0; i < mesh_->n_numbers_; ++i) {
    auto theta = i * d_theta;
    auto x     = mesh_->chord_ / 2.0 + mesh_->radius_ * std::cos(theta);
    auto y     = -mesh_->radius_ * std::sin(theta);
    mesh_->h_points_(mesh_->m_numbers_ - 1, i, 0) = x;
    mesh_->h_points_(mesh_->m_numbers_ - 1, i, 1) = y;
  }
  for (int i = 1; i < mesh_->m_numbers_; ++i) {
    auto surface_subview =
        Kokkos::subview(mesh_->h_points_, 0, Kokkos::ALL, Kokkos::ALL);
    auto far_subview = Kokkos::subview(mesh_->h_points_, mesh_->m_numbers_ - 1,
                                       Kokkos::ALL, Kokkos::ALL);
    for (int j = 0; j < mesh_->n_numbers_; ++j) {
      mesh_->h_points_(i, j, 0) =
          surface_subview(j, 0) +
          i * (far_subview(j, 0) - surface_subview(j, 0)) /
              (mesh_->m_numbers_ - 1.0);
      mesh_->h_points_(i, j, 1) =
          surface_subview(j, 1) +
          i * (far_subview(j, 1) - surface_subview(j, 1)) /
              (mesh_->m_numbers_ - 1.0);
    }
  }
  Kokkos::deep_copy(mesh_->points_, mesh_->h_points_);
}

std::shared_ptr<MeshSolverBase> MeshSolverO::create(
    const MeshParams::MeshParamsPtr& mesh_ptr) {
  auto solver = std::make_shared<MeshSolverO>(mesh_ptr);
  return solver;
}

MeshSolverFactory::MeshSolverFactory() {
  registerSolver("O", &MeshSolverO::create);
}

MeshSolverFactory::~MeshSolverFactory() { solver_map_.clear(); }

MeshSolverFactory& MeshSolverFactory::Instance() {
  static MeshSolverFactory instance;
  return instance;
}

void MeshSolverFactory::registerSolver(const std::string& key,
                                       SolverCreateFn cfn) {
  solver_map_[key] = cfn;
}

std::shared_ptr<MeshSolverBase> MeshSolverFactory::create(
    const MeshParams::MeshParamsPtr& mesh_params) {
  auto it = solver_map_.find(mesh_params->type_);
  if (it != solver_map_.end()) {
    auto mesh_solver = it->second(mesh_params);
    mesh_solver->init();
    return mesh_solver;
  }
  fmt::print("#########ERROR#########\n");
  fmt::print("This Code do not implement the {} mesh solver\n",
             mesh_params->type_);
  return nullptr;
}

}  // namespace cfd_kokkos::mesh
