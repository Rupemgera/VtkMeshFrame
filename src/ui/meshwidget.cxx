#include "meshwidget.h"
#include "ui_meshwidget.h"
#include <fstream>
#include <iostream>
#include <qcolordialog.h>
#include <qmessagebox.h>
#include <sstream>

void
MeshWidget::initialization(QWidget* parent,
                           const json11::Json& json,
                           std::string dir_path)
{
  ui->setupUi(this);
  _viewer = std::make_unique<viewtools::ViewManager>(ui->viewerWidget);
  _shell = std::make_unique<MeshShell>(_viewer.get());

  if (dir_path != "") {
    dir_path += "/";
  }
  // mesh
  auto mesh_file = json["mesh"].string_value();
  if (mesh_file != "") {
    ui->pushButton_read->setDisabled(true);

    _shell->readMesh(dir_path + mesh_file);
    _shell->drawMesh(getRenderStyle());
    updateMeshOpacity();
    updateMeshInfo();
  }

  addSlot();
};

MeshWidget::MeshWidget(QWidget* parent)
  : QWidget(parent)
  , ui(new Ui::MeshWidget)
{
  initialization(parent, json11::Json{});
}

MeshWidget::~MeshWidget()
{
  delete ui;
}

void
MeshWidget::addSlot()
{}

void
MeshWidget::updateMeshInfo()
{
  std::string info = _shell->mesh_info();
  ui->textBrowser_mesh->append(info.c_str());
  // std::cout << info << std::endl;
}

int
MeshWidget::getRenderStyle()
{
  int nRenderStyle = 0;
  if (ui->checkBox_edge->isChecked() == true) {
    nRenderStyle += 1;
  }
  if (ui->checkBox_face->isChecked() == true) {
    nRenderStyle += 2;
  }
  return nRenderStyle;
}

std::string
MeshWidget::filenameFromDialog(const char* dialog_name, const char* filter)
{
  QString qfname =
    QFileDialog::getOpenFileName(0, dialog_name, _directory_path, filter);

  if (qfname == "") {
    return "";
  }

  /********** save last opened path *************/

  int i = qfname.lastIndexOf('/');
  _directory_path = qfname.left(i);

  return qfname.toStdString();
}

void
MeshWidget::messageBox(const char* info)
{
  QMessageBox::information(NULL, "Warning", info);
}

void
MeshWidget::readMesh()
{
  std::string filename = filenameFromDialog(
    "Open Mesh File", "OVM files(*.ovm);;Abaqus inp files(*.inp)");

  ui->pushButton_read->setDisabled(true);

  /********** Mesh *************/

  _shell->readMesh(filename);
  _shell->drawMesh(getRenderStyle());
  updateMeshOpacity();
  updateMeshInfo();
}

void
MeshWidget::updateMeshRenderStyle()
{
  _shell->updateMeshRenderStyle(getRenderStyle());
}

void
MeshWidget::updateMeshOpacity()
{
  double opacity = ui->doubleSpinBox_opacity->value();
  _shell->updateFaceOpacity(opacity);
}

void
MeshWidget::on_pushButton_refresh_listWidget_clicked()
{
  auto additem3d = [&](viewtools::Actor3dControler* ac) {
    QListWidgetItem* new_item = new QListWidgetItem(ac->getName().c_str());
    /*new_item->setFlags(Qt::ItemIsSelectable | Qt::ItemIsUserCheckable |
                       Qt::ItemIsEnabled);*/
    new_item->setFlags(new_item->flags() | Qt::ItemIsUserCheckable);

    if (ac->getVisibility())
      new_item->setCheckState(Qt::Checked);
    else
      new_item->setCheckState(Qt::Unchecked);
    ui->listWidget_actor3d->addItem(new_item);
  };

  auto additem2d = [&](viewtools::Actor2dControler* ac) {
    QListWidgetItem* new_item = new QListWidgetItem(ac->getName().c_str());
    /*new_item->setFlags(Qt::ItemIsSelectable | Qt::ItemIsUserCheckable |
                       Qt::ItemIsEnabled);*/
    new_item->setFlags(new_item->flags() | Qt::ItemIsUserCheckable);

    if (ac->getVisibility())
      new_item->setCheckState(Qt::Checked);
    else
      new_item->setCheckState(Qt::Unchecked);
    ui->listWidget_actor2d->addItem(new_item);
  };

  ui->listWidget_actor2d->clear();
  ui->listWidget_actor3d->clear();

  for (auto& ac : _viewer->getTable3d()) {
    // mesh render is not added
    if (ac.second->getName().find("mesh") == ac.second->getName().npos)
      additem3d(ac.second.get());
  }

  for (auto& ac : _viewer->getTable2d()) {
    additem2d(ac.second.get());
  }
}

void
MeshWidget::on_listWidget_actor3d_itemChanged(QListWidgetItem* item)
{
  std::string name = item->text().toStdString();
  if (item->checkState() == Qt::Checked) {
    _viewer->setVisibility(name, true);
  } else {
    _viewer->setVisibility(name, false);
  }
  _viewer->refresh();
}

void
MeshWidget::on_listWidget_actor3d_itemClicked(QListWidgetItem* item)
{
  std::string name = item->text().toStdString();
  auto status = _viewer->getStatus(name);
  ui->actor_opacity->setValue(std::get<0>(status));
  ui->actor_size->setValue(std::get<1>(status));

  auto c = _viewer->getColor(name);
  QColor color(std::get<0>(c), std::get<1>(c), std::get<2>(c));
  ui->color_frame->setPalette(QPalette(color));
  // std::cout << name << " clicked" << std::endl;
}

void
MeshWidget::on_listWidget_actor2d_itemChanged(QListWidgetItem* item)
{
  std::string name = item->text().toStdString();
  if (item->checkState() == Qt::Checked) {
    _viewer->setVisibility(name, true);
  } else {
    _viewer->setVisibility(name, false);
  }
  _viewer->refresh();
}

void
MeshWidget::on_listWidget_actor2d_itemClicked(QListWidgetItem* item)
{}

void
MeshWidget::on_pushButton_actor_color_clicked()
{
  QColor color = QColorDialog::getColor();
  if (color.isValid()) {
    // set actor's color
    double red = color.red();
    double green = color.green();
    double blue = color.blue();
    double c[] = { red / 255.0, green / 255.0, blue / 255.0 };
    auto item = ui->listWidget_actor3d->currentItem();
    std::string name = item->text().toStdString();
    _viewer->setColor(name, c);

    // gui changes
    ui->color_frame->setPalette(QPalette(color));
  }
}

void
MeshWidget::on_actor_refresh_clicked()
{
  double size = ui->actor_size->value();
  double opacity = ui->actor_opacity->value();
  auto item = ui->listWidget_actor3d->currentItem();
  std::string name = item->text().toStdString();
  _viewer->setSize(name, size);
  _viewer->setOpacity(name, opacity);
}

void
MeshWidget::on_checkBox_vertex_id_toggled()
{
  //_shell->showIds();
}
