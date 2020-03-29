#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget* parent)
  : QMainWindow(parent)
  , ui(new Ui::MainWindow)
{
  ui->setupUi(this);

  addSlot();
}

MainWindow::~MainWindow()
{
  delete ui;
}

void
MainWindow::insert_tabWidget(MeshWidget* widget, QString widget_name)
{
  ui->tabWidget->addTab(widget, "Mesh Shell");
  ui->tabWidget->setCurrentWidget(widget);
  ui->tabWidget->show();
}

// void MainWindow::browserPrint(QString text)
// {
// 	//this->ui->textBrowser->append(text + "\n");

// 	//ConsoleWidget* con =
// 	_console_widget->print(text);
// }

void
MainWindow::addSlot()
{
  connect(this->ui->tabWidget,
          SIGNAL(tabCloseRequested(int)),
          this,
          SLOT(closeMeshShell(int)));
}

void
MainWindow::closeMeshShell(int id)
{
  auto shell = ui->tabWidget->widget(id);
  ui->tabWidget->removeTab(id);
  delete shell;
}

void
MainWindow::on_actionNewWidget_triggerd()
{
  auto mesh_widget = new MeshWidget(this);
  insert_tabWidget(mesh_widget, "Empty MeshShell");
}
