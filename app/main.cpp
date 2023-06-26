#include <iostream>
#include <Kokkos_Core.hpp>
#include "util/parameters.hpp"
#include "mesh/solver.hpp"

int main(int argc, char *argv[]) {
  std::string config_file;
  if (argc < 1) {
    const std::string path   = argv[0];
    std::size_t last_sep_pos = path.find_last_of("/\\");
    std::string directory    = path.substr(0, last_sep_pos);
    directory.append("/config.json");
    config_file = directory;
  } else {
    config_file = argv[1];
  }
  Kokkos::initialize(argc, argv);
  {
    auto parameters  = std::make_shared<cfd_kokkos::Params>(config_file);
    auto mesh_solver = cfd_kokkos::mesh::MeshSolverFactory::Instance().create(
        parameters->mesh_params_ptr_);
    mesh_solver->solve();
    mesh_solver->save_mesh();
  }
  Kokkos::finalize();
  return 0;
}  // end main