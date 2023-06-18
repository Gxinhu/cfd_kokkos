# Two alternatives:
# 1. If KOKKOS_PROJ_TMPL_BUILD_KOKKOS is ON, we download kokkos sources and build them using FetchContent (which actually uses add_subdirectory)
# 2. If KOKKOS_PROJ_TMPL_BUILD_KOKKOS is OFF (default), we don't build kokkos, but use find_package for setup (you must have kokkos already installed)

# NOTE about required C++ standard
# we better chose to set the minimum C++ standard level if not already done:
# - when building kokkos <  4.0.00, it defaults to c++-14
# - when building kokkos >= 4.0.00, it defaults to c++-17
# - when using installed kokkos, we set C++ standard according to kokkos version

#
# Do we want to build kokkos (https://github.com/kokkos/kokkos) ?
#

#
# Option to use git (instead of tarball release) for downloading kokkos
#

option(KOKKOS_PROJ_TMPL_USE_GIT "Turn ON if you want to use git to download Kokkos sources (default: OFF)" ON)


# find_package(Git REQUIRED)
include(FetchContent)

if(KOKKOS_PROJ_TMPL_USE_GIT)
  FetchContent_Declare(fmt_external
    GIT_REPOSITORY https://gitee.com/chooosky/fmt.git
    GIT_TAG 9.1.0
  )
else()
  FetchContent_Declare(fmt_external
    URL https://github.com/fmtlib/fmt/archive/refs/tags/fmt-10.0.0.zip
  )
endif()

# Import kokkos targets (download, and call add_subdirectory)
FetchContent_MakeAvailable(fmt_external)
