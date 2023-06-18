#include "mesh/solver.hpp"

namespace cfd_kokkos::mesh {

MeshSolverBase::MeshSolverBase(const MeshParams::MeshPtr& mesh_ptr) {
  mesh_ = cfd_kokkos::mesh::MeshAirfoilFactory::create(mesh_ptr);
}

MeshSolverO::MeshSolverO(const MeshParams::MeshPtr& mesh_ptr)
    : MeshSolverBase(mesh_ptr) {}

void MeshSolverO::solve() { fmt::print("Hello World\n"); }

void MeshSolverO::init() { fmt::print("Hello World\n"); }

std::shared_ptr<MeshSolverBase> MeshSolverO::create(
    const MeshParams::MeshPtr& mesh_ptr) {
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
    const MeshParams::MeshPtr& mesh_params) {
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
