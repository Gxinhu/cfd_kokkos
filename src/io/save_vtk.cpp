#include "io/save_vtk.hpp"
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkXMLDataSetWriter.h>
#include <cstddef>
#include "mesh/mesh_shared.hpp"
#include "util/cfd_shared.hpp"

namespace cfd_kokkos::io {
void save_vtk_2d(const mesh::MeshMatrix::HostMirror &h_points,
                 const std::string &save_name) {
  const auto nx = h_points.extent(0);
  const auto ny = h_points.extent(1);

  vtkSmartPointer<vtkStructuredGrid> structured_grid =
      vtkSmartPointer<vtkStructuredGrid>::New();
  structured_grid->SetDimensions(nx, ny, 1);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(nx * ny);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      auto x = h_points(i, j, 0);
      auto y = h_points(i, j, 1);
      auto z = 0.0;
      points->SetPoint(i * ny + j, x, y, z);
    }
  }
  structured_grid->SetPoints(points);
  vtkSmartPointer<vtkXMLDataSetWriter> writer =
      vtkSmartPointer<vtkXMLDataSetWriter>::New();
  writer->SetInputData(structured_grid);
  writer->SetFileName(save_name.c_str());
  writer->Write();
}

}  // namespace cfd_kokkos::io