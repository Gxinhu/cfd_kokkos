#pragma once
#include <cstddef>
#include <memory>
#include "cfd_shared.hpp"
#include "mesh/mesh_shared.hpp"
#include <fstream>
namespace cfd_kokkos {

struct MeshParams {
  using MeshParamsPtr = std::shared_ptr<MeshParams>;
  precision chord_{};
  precision thickness_{};
  precision max_camber_{};
  precision max_location_{};
  size_t n_numbers_{};
  size_t m_numbers_{};
  size_t display_iter_{};
  size_t max_iter_{};
  precision radius_{};
  precision error_{};
  std::string type_{};
  bool is_sys_ = false;
  void setup(json& data) {
    chord_        = data["chord"].get<precision>();
    thickness_    = data["thickness"].get<precision>();
    n_numbers_    = data["n_numbers"].get<int>();
    m_numbers_    = data["m_numbers"].get<int>();
    max_camber_   = data["max_camber"].get<int>();
    max_location_ = data["max_location"].get<int>();
    display_iter_ = data["display_iter"].get<int>();
    max_iter_ = data["max_iter"].get<int>();
    is_sys_       = data["is_sys"].get<bool>();
    radius_       = data["radius"].get<precision>();
    error_       = data["error"].get<precision>();
    type_         = data["type"].get<std::string>();
  }
};
class Params {
 public:
  explicit Params(const std::string& config_file) {
    std::ifstream f;
    std::cout << "Load Config File " << config_file << '\n';
    f.open(config_file);
    if (!f) {
      throw std::runtime_error("Could not Open file");
    }
    auto data = json::parse(f);
    mesh_params_ptr_->setup(data["mesh"]);
    data.clear();
    f.close();
  };

  MeshParams::MeshParamsPtr mesh_params_ptr_ = std::make_shared<MeshParams>();
};
}  // namespace cfd_kokkos