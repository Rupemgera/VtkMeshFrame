#pragma once
#include "mesh/mesh_wrapper.h"
#include "view_tools/vtk_wrapper.h"

/*********** defines begin **************/

using VtkWrapper = viewtools::VtkWrapper;

struct MeshData
{
  std::vector<Eigen::Vector3d> points;
  std::vector<Eigen::Matrix<long long, 3, 1>> faces;
  std::vector<Eigen::Matrix<long long, 3, 1>> boundary_faces;
  void getFaceData(const std::vector<int>& ids,
                   std::vector<Eigen::Matrix<long long, 3, 1>>& boundary_faces)
  {
    boundary_faces.clear();
    for (int i : ids) {
      boundary_faces.push_back(faces[i]);
    }
  }
};

class MeshShell
{
protected:
  VtkWrapper* _viewer;

  std::unique_ptr<MeshWrapper> mesh_wrapper = nullptr;

  std::string mesh_name = "mesh";

public:
  MeshShell(VtkWrapper* viewer);
  ~MeshShell();

  void drawMesh(int nRenderStyle = 3);

  void readMesh(std::string filename);

  std::string mesh_info() const;

  void updateMeshRenderStyle(int nRenderStyle);

  /**
   *geometry = 1 : normal mesh
   *geometry = 2 : shrinked mesh
   */
  void updateFaceOpacity(double opacity);

  /********** data ************/
  MeshData mesh_data;
};