#pragma once

#include "json11/json11.hpp"
#include "view_tools/view_manager.h"
#include "mesh_shell.h"

#include "QtWidgets/qlistwidget.h"
#include <QFile>
#include <QFileDialog>
#include <QWidget>
#include <qstandardpaths.h>

#include <cstdarg>
#include <initializer_list>
#include <map>
#include <memory>

using ActorMap = std::map<std::string, std::string>;

namespace Ui {
class MeshWidget;
}

class MeshWidget : public QWidget
{
  Q_OBJECT
private:
  void initialization(QWidget* parent,
                      const json11::Json& json,
                      std::string dir_path = "");

public:
  MeshWidget(QWidget* parent);
  ~MeshWidget();

private:
  Ui::MeshWidget* ui;

  void addSlot();

  std::unique_ptr<viewtools::ViewManager> _viewer;
  std::unique_ptr<MeshShell> _shell;
  ActorMap active_actors;

  bool _mesh_loaded = false;

  /* functions */

  void updateMeshInfo();

  int getRenderStyle();

  std::string filenameFromDialog(const char* dialog_name, const char* filter);

  void messageBox(const char* info);

  /* properties */
#ifdef __linux
  QString _directory_path = "~/";
#endif
#ifdef WIN32
  QString _directory_path =
    QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation) +
    "/Models";
#endif

private slots:
  void readMesh();

  void updateMeshRenderStyle();

  void updateMeshOpacity();


  /******************** listWidget ******************************/

  void on_pushButton_refresh_listWidget_clicked();

  void on_listWidget_actor3d_itemChanged(QListWidgetItem* item);

  void on_listWidget_actor3d_itemClicked(QListWidgetItem* item);

  void on_listWidget_actor2d_itemChanged(QListWidgetItem* item);

  void on_listWidget_actor2d_itemClicked(QListWidgetItem* item);

  void on_pushButton_actor_color_clicked();

  void on_actor_refresh_clicked();

  void on_checkBox_vertex_id_toggled();


 };