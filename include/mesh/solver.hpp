#pragma once

#include <memory>
#include <shared_mutex>
#include "airfoil.hpp"
#include "fmt/core.h"
#include "mesh/mesh_shared.hpp"
#include "util/parameters.hpp"
namespace cfd_kokkos::mesh {

class MeshSolverBase {
 public:
  explicit MeshSolverBase(const MeshParams::MeshParamsPtr& mesh_ptr);
  virtual void solve() = 0;
  virtual void init()  = 0;
  std::shared_ptr<NACA4DigitBase> mesh_;
  void save_mesh() const;
};

class MeshSolverO : public MeshSolverBase {
 public:
  explicit MeshSolverO(const MeshParams::MeshParamsPtr& mesh_ptr);
  void solve() override;
  void init() override;
  static std::shared_ptr<MeshSolverBase> create(
      const MeshParams::MeshParamsPtr& mesh_ptr);
};

class MeshSolverC {};

class MeshSolverPossion {};

class MeshSolverFactory {
 private:
  MeshSolverFactory();
  MeshSolverFactory& operator=(const MeshSolverFactory&) {
    return *this;
  }  // non-copyable

  using SolverCreateFn = std::shared_ptr<MeshSolverBase> (*)(
      const MeshParams::MeshParamsPtr& mesh_params);
  using SolverCreateMap = std::map<std::string, SolverCreateFn>;
  SolverCreateMap solver_map_;

 public:
  ~MeshSolverFactory();
  MeshSolverFactory(const MeshSolverFactory&) = delete;
  static MeshSolverFactory& Instance();
  void registerSolver(const std::string& key, SolverCreateFn cfn);
  std::shared_ptr<MeshSolverBase> create(
      const MeshParams::MeshParamsPtr& mesh_params);
};

}  // namespace cfd_kokkos::mesh