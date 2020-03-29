#include "mesh_shell.h"

MeshShell::MeshShell(VtkWrapper* viewer)
  : _viewer(viewer)
{
  mesh_wrapper = std::make_unique<MeshWrapper>();
}

/**
 *viewer should be released by caller
 */
MeshShell::~MeshShell() {}

void
MeshShell::drawMesh(int nRenderStyle)
{
  if (!_viewer->exist(mesh_name)) {
    std::vector<int> boundary_face_id_list;

    mesh_wrapper->get_boundary_faces(boundary_face_id_list);
    mesh_wrapper->get_face_data<3>(mesh_data.points, mesh_data.faces);
    mesh_data.getFaceData(boundary_face_id_list, mesh_data.boundary_faces);

    mesh_name = "mesh";

    _viewer->drawTetMesh(mesh_name, mesh_data.points, mesh_data.boundary_faces);

    // first render, reset camera
    _viewer->resetCamera();
  } else {
    _viewer->setVisibility(mesh_name, true);
  }
  _viewer->setRenderStyle(mesh_name, nRenderStyle);
  _viewer->refresh();
}

void
MeshShell::readMesh(std::string filename)
{
  mesh_wrapper->read_mesh(filename);
}

std::string
MeshShell::mesh_info() const
{
  return mesh_wrapper->mesh_info();
}

void
MeshShell::updateMeshRenderStyle(int nRenderStyle)
{
  _viewer->setRenderStyle(mesh_name, nRenderStyle);
  _viewer->refresh();
}

void
MeshShell::updateFaceOpacity(double opacity)
{
  _viewer->setOpacity(mesh_name, opacity);
}
