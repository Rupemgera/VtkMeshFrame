#pragma once
/**
 * @file view_manager.h
 * @author Liao Yizhou (3130001008@zju.edu.cn)
 * @brief manager actors and transport data to widgets
 * @version 0.1
 * @date 2020-02-20
 *
 * @copyright Copyright (c) 2020
 *
 */

#include "actor_controler.h"
#include "vtk_wrapper.h"

namespace viewtools {

class ViewManager : public VtkWrapper {
protected:
  /** Axes
   * add a Orthogonal coordinate to the bottom left of the screen
   */
  void addAxes();

  /************ vtk renderer ************************/

  vtkSmartPointer<vtkRenderer> _renderer;
  // renderWindow
  vtkSmartPointer<vtkGenericOpenGLRenderWindow> _renderWindow;
  vtkSmartPointer<vtkOrientationMarkerWidget> _markerWidget;
  QVTKOpenGLWidget *_vtkWidget;

  /******************** source management ***************/

  // when remove an ActorControler from the table, the actor should also be
  // removed from the renderer
  // for Actor2dPtr's custom deleter call
  // removeActor, which call _renderer.removeActor, so actor_tables must be
  // declared after _renderer. Thus tables' destructors will be called before
  // _renderer's destructor.
  std::function<void(Actor2dControler *)> actor2d_deleter =
      [this](Actor2dControler *ptr) {
        this->removeActor(ptr->getProp());
        delete ptr;
      };
  std::function<void(Actor3dControler *)> actor3d_deleter =
      [this](Actor3dControler *ptr) {
        this->removeActor(ptr->getProp());
        delete ptr;
      };

  std::map<std::string, Actor2dPtr> _table2d_;

  std::map<std::string, Actor3dPtr> _table3d_;

public:
  ViewManager(QVTKOpenGLWidget *widget);

  const std::map<std::string, Actor2dPtr> &getTable2d() const {
    return _table2d_;
  };

  const std::map<std::string, Actor3dPtr> &getTable3d() const {
    return _table3d_;
  }

  void refresh() override;

  /* reset camera */
  void resetCamera() override;

protected:
  /************** source management ********************/

  /**
   * @brief render the actor
   *
   * @param actor
   * @return int numbers of actors in the renderer
   */
  int addActor(vtkProp *actor);

  /**
   * @brief remove the actor
   *
   * @param actor
   * @return int numbers of actors in the renderer
   */
  int removeActor(vtkProp *actor);

  void insert(Actor2dControler *ac) override;
  void insert(Actor3dControler *ac) override;
  void remove(std::string name) override;
  std::tuple<bool, Actor2dControler *> findActor2d(std::string name) override;
  std::tuple<bool, Actor3dControler *> findActor3d(std::string name) override;
};
} // namespace viewtools