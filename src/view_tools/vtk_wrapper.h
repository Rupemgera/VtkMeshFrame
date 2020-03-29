/**
 * @file vtk_wrapper.h
 * @author Liao Yizhou (3130001008@zju.edu.cn)
 * @brief VtkWrapper
 * @version 0.1
 * @date 2020-01-03
 *
 * @copyright Copyright (c) 2020
 *
 */
#pragma once
#include "actor_controler.h"
#include "view_defs.h"
#include "vtk_templates.hpp"
#include <functional>
#include <map>
#include <memory>
#include <tuple>
#include <vector>

namespace viewtools {
using Actor2dPtr =
  std::unique_ptr<Actor2dControler, std::function<void(Actor2dControler*)>>;
using Actor3dPtr =
  std::unique_ptr<Actor3dControler, std::function<void(Actor3dControler*)>>;

/******************* defines ***************************/
class VtkWrapper
{
protected:
  /*********** process data begin **************************/

  template<typename T, int cell_n>
  vtkSmartPointer<vtkPolyData> processPolyData(
    const std::vector<Eigen::Vector3d>& points,
    const std::vector<Eigen::Matrix<T, cell_n, 1>>& polys);

  vtkSmartPointer<vtkPolyData> processPolyData(
    const RowMatrix<double>& points,
    const RowMatrix<IndexType>& polys);

  vtkSmartPointer<vtkPolyData> processPointData(
    const RowMatrix<double>& points);
  /**
   *@param name is the key for seaching
   */
  template<int cell_n>
  bool drawMesh(std::string name,
                const std::vector<Eigen::Vector3d>& points,
                const std::vector<VertexList<cell_n>>& faces);

  /**
   * @brief convert points in Eigen matrix form to VTK array
   *
   * @param points coordinates of points, (x,y,z) stored in row. total size is
   * n*3.
   * @return vtkSmartPointer<vtkPoints>
   */
  vtkSmartPointer<vtkPoints> convert_points(const RowMatrix<double>& points);

  vtkSmartPointer<vtkPoints> convert_points(
    const std::vector<Eigen::Vector3d>& points);

  vtkSmartPointer<vtkDoubleArray> convert_arrary(
    const RowMatrix<double>& matrix);

  vtkSmartPointer<vtkDoubleArray> convert_arrary(
    const std::vector<double>& scalars);

  vtkSmartPointer<vtkDoubleArray> convert_vectors(
    const std::vector<Eigen::Vector3d>& vecs);

  vtkSmartPointer<vtkDoubleArray> convert_vectors(
    const RowMatrix<double>& vecs);

  std::vector<vtkSmartPointer<vtkDoubleArray>> convert_frames(
    const std::vector<Eigen::Matrix<double, 3, 3>>& vecs);

  /***************** source management *********************/

  // virtual int addActor(vtkProp *actor) = 0;

  // virtual int removeActor(vtkProp *actor) = 0;

  virtual void insert(Actor2dControler* ac) = 0;
  virtual void insert(Actor3dControler* ac) = 0;
  virtual void remove(std::string name) = 0;
  virtual std::tuple<bool, Actor2dControler*> findActor2d(std::string name) = 0;
  virtual std::tuple<bool, Actor3dControler*> findActor3d(std::string name) = 0;

  /***************** color ************************/

  /* color arrays */
  std::vector<Color> _preset_colors;

  void init_preset_colors();

public:
  VtkWrapper() { init_preset_colors(); }

  virtual ~VtkWrapper() = default;

  // virtual const std::map<std::string, Actor3dPtr> &getTable3d() const = 0;
  // virtual const std::map<std::string, Actor2dPtr> &getTable2d() const = 0;

  virtual void refresh() = 0;
  virtual void resetCamera() = 0;

  // test interface
  void testRenderFunction();

  bool drawArrows(std::string name,
                  const RowMatrix<double>& points,
                  const RowMatrix<double>& vectors);

  bool drawArrows(std::string name,
                  const std::vector<Eigen::Vector3d>& points,
                  const std::vector<Eigen::Vector3d>& vectors);

  bool drawFaces(std::string name,
                 const RowMatrix<double>& points,
                 const RowMatrix<IndexType>& faces);

  bool drawFrames(std::string name,
                  const RowMatrix<double>& points,
                  const std::vector<Eigen::Matrix<double, 3, 3>>& frames,
                  double scale_factor = 0.8,
                  double line_width = 1.0);

  bool drawGraph(
    std::string name,
    const std::vector<Eigen::Vector3d>& points,
    const std::vector<Eigen::Matrix<IndexType, 2, 1>>& vertices_pairs,
    std::vector<double>* scalars = nullptr,
    double line_width = 2.0);

  bool drawGraph(std::string name,
                 const RowMatrix<double>& points,
                 const RowMatrix<IndexType>& edge_pairs,
                 const RowMatrix<double>* scalars = nullptr,
                 double line_width = 2.0);

  bool showIds(std::string name);

  bool drawIsosurface(std::string name,
                      const RowMatrix<double>& points,
                      const RowMatrix<IndexType>& faces,
                      const RowMatrix<double>& scalars,
                      double value);

  /**
   * @brief draw points with labels
   *
   * @param name
   * @param points matrix of n*3, presenting n points
   * @param labels matrix of n*1, presenting points' labels
   * @return true
   * @return false
   */
  bool drawLabels(std::string name,
                  const RowMatrix<double>& points,
                  const RowMatrix<int>& labels);

  bool drawLabels(std::string name,
                  const RowMatrix<double>& points,
                  const std::vector<int>* labels = nullptr);

  bool drawLabelsFloat(std::string name,
                       const RowMatrix<double>& points,
                       const std::vector<double>* labels = nullptr);

  /**
   * @brief draw segment lines
   */
  bool drawLines(std::string name,
                 const std::vector<std::vector<Eigen::Vector3d>>& points,
                 bool is_loop = false);

  bool drawMolecule(std::string name,
                    const RowMatrix<double>& points,
                    const RowMatrix<IndexType>& bonds,
                    const RowMatrix<IndexType>* scalars = nullptr);

  bool drawMoledule(std::string name,
                    const RowMatrix<double>& points,
                    std::vector<unsigned short>* label = nullptr);

  bool drawPoints(std::string name,
                  const std::vector<Eigen::Vector3d>& points,
                  double point_size,
                  std::vector<unsigned short>* label = nullptr);

  bool drawPoints(std::string name,
                  const RowMatrix<double>& points,
                  double point_size = 2.0,
                  std::vector<unsigned short>* label = nullptr);

  bool drawSegments(
    std::string name,
    const std::vector<Eigen::Vector3d>& points,
    const std::vector<Eigen::Matrix<IndexType, 2, 1>>& vertices_pairs,
    double line_width = 2.0);

  bool drawTetMesh(std::string name,
                   const std::vector<Eigen::Vector3d>& points,
                   const std::vector<Eigen::Matrix<IndexType, 3, 1>>& faces);

  /* vector field */

  bool drawVector(std::string name,
                  const std::vector<Eigen::Vector3d>& points,
                  const std::vector<Eigen::Vector3d>& vectors,
                  double scale_factor = 0.8,
                  double line_width = 1.0);

  bool drawVector(std::string name,
                  const RowMatrix<double>& points,
                  const RowMatrix<double>& vectors,
                  double scale_factor = 1.0,
                  double line_width = 1.0);

  /**
   * @brief draw lines that represent vectors. the lines are colored according
   * to the scalars of the vectors.
   *
   * @param name
   * @param points
   * @param vectors
   * @param scalars
   * @param scale_factor the length of the lines. if equals 1.0, then the length
   * will be as vectors' norm.
   * @param line_width
   * @return true
   * @return false
   */
  bool drawVectorWithScalars(std::string name,
                             const std::vector<Eigen::Vector3d>& points,
                             const std::vector<Eigen::Vector3d>& vectors,
                             std::vector<double>& scalars,
                             double scale_factor = 0.8,
                             double line_width = 1.0);

  /**
   *@param opacity from 0.0 to 1.0
   */
  void setOpacity(std::string name, double opacity);

  void setVisibility(std::string name, bool flag);

  void setColor(std::string name, double* color);

  std::tuple<int, int, int> getColor(std::string name);

  void setSize(std::string name, double size);

  double getSize(std::string name);

  // 1.opacity 2.size
  std::tuple<double, double> getStatus(std::string name);
  /**
   *@param render_style 1: only edges 2: only faces 3:edges and faces
   */
  void setRenderStyle(std::string name, int render_style);

  /* Vertex Scalar */
  void setVertexScalars(std::string name,
                        std::vector<double>& scalars,
                        double lower_bound,
                        double upper_bound);

  bool exist(std::string name);
};
} // namespace viewtools

#include "vtk_wrapper_impl.hpp"