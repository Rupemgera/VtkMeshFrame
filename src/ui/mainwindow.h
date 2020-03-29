#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "meshwidget.h"

#include <memory>

#include <QFile>
#include <QFileDialog>
#include <QMainWindow>
#include <QWidget>
#include <qdir.h>
#include <qstandardpaths.h>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  explicit MainWindow(QWidget* parent = nullptr);
  ~MainWindow();

private:
  /*
  Properties
  */

  /*UI*/

  Ui::MainWindow* ui;

  /*
  Widgets
  */

  void insert_tabWidget(MeshWidget* widget, QString widget_name);

  // MeshWidget *_mesh_widget;

  // vtkSmartPointer<vtkInteractorStyle> _style;

  /*
  Functions
  */

  /*widget functions*/

  // void browserPrint(QString text);

  void addSlot();

  /*slots functions*/
private slots:

  void closeMeshShell(int id);

  void on_actionNewWidget_triggerd();
};

#endif // MAINWINDOW_H
