#include "view_manager.h"

namespace viewtools {

ViewManager::ViewManager(QVTKOpenGLWidget *widget) {
  // Create the usual rendering stuff

  _renderer = vtkRenderer::New();

  _renderWindow = vtkGenericOpenGLRenderWindow::New();
  _renderWindow->AddRenderer(_renderer);

  _vtkWidget = widget;
#if VTK_MAJOR_VERSION == 8 && VTK_MINOR_VERSION >= 9
  _vtkWidget->setRenderWindow(_renderWindow);
#else
  _vtkWidget->SetRenderWindow(_renderWindow);
#endif
  //_renderWindow = (vtkGenericOpenGLRenderWindow)_vtkWidget->renderWindow();
  // interaction style
  vtkNew<vtkInteractorStyleSwitch> style_switch;
  vtkNew<vtkInteractorStyleTrackballCamera> style_trackball_camera;
  _renderWindow->GetInteractor()->SetInteractorStyle(style_trackball_camera);

  addAxes();
}

void ViewManager::refresh() { _renderWindow->Render(); }

void ViewManager::resetCamera() { _renderer->ResetCamera(); }

void ViewManager::addAxes() {
  vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
  _markerWidget = vtkOrientationMarkerWidget::New();
  //_markerWidget->SetOutlineColor(0.9300, 0.5700, 0.1300);
  _markerWidget->SetOrientationMarker(axes);
  _markerWidget->SetInteractor(_renderWindow->GetInteractor());
  _markerWidget->SetViewport(0.0, 0.0, 0.2, 0.2);
  _markerWidget->SetEnabled(1);
  _markerWidget->InteractiveOff();
}

int ViewManager::addActor(vtkProp *actor) {
  // Add Actor to renderer
  _renderer->AddActor(actor);
  //_renderer->ResetCamera();
  _renderWindow->Render();
  //_vtkWidget->show();

  return _renderer->GetActors()->GetNumberOfItems();
}

int ViewManager::removeActor(vtkProp *actor) {
  _renderer->RemoveActor(actor);
  _renderWindow->Render();
  return _renderer->GetActors()->GetNumberOfItems();
}

void ViewManager::insert(Actor3dControler *ac) {
  // render actor
  addActor(ac->getProp());
  std::string name = ac->getName();
  _table3d_[name] = Actor3dPtr(ac, actor3d_deleter);
}

void ViewManager::insert(Actor2dControler *ac) {
  addActor(ac->getProp());
  std::string name = ac->getName();
  _table2d_[name] = Actor2dPtr(ac, actor2d_deleter);
}

void ViewManager::remove(std::string name) {
  auto target = _table3d_.find(name);
  if (target != _table3d_.end()) {
    _table3d_.erase(target);
  } else {
    auto t2 = _table2d_.find(name);
    if (t2 != _table2d_.end()) {
      _table2d_.erase(name);
    }
  }
}
std::tuple<bool, Actor2dControler *>
ViewManager::findActor2d(std::string name) {
  auto target = _table2d_.find(name);
  if (target != _table2d_.end()) {
    return std::make_tuple(true, target->second.get());
  } else {
    return std::make_tuple(false, nullptr);
  }
}
std::tuple<bool, Actor3dControler *>
ViewManager::findActor3d(std::string name) {
  auto target = _table3d_.find(name);
  if (target != _table3d_.end()) {
    return std::make_tuple(true, target->second.get());
  } else {
    return std::make_tuple(false, nullptr);
  }
}
} // namespace viewtools