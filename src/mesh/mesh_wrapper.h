#pragma once
/**
 * @file mesh_wrapper.h
 * @author Liao Yizhou (3130001008@zju.edu.cn)
 * @brief class MeshWrapper
 * @version 0.1
 * @date 2019-09-08
 *
 * @copyright Copyright (c) 2019
 *
 */

#include "../common.h"
#include "meshDefs.h"

/*
data struct defines
*/

// TF is a struct to store 3 vertices of a halfface
// the order of vertices is ignored,
// which means hf1 = 1-2-3 equals hf2 = 2-3-1
// while hf1 = 1-2-3 does not equals hf2 = 3-2-1
struct TF
{
  int first;
  int second;
  int third;

  TF(int x, int y, int z)
  {
    std::vector<int> ve;
    ve.push_back(x), ve.push_back(y), ve.push_back(z);
    std::sort(ve.begin(), ve.end());
    first = ve[0], second = ve[1], third = ve[2];
  }
  bool operator<(const TF& tf2) const
  {
    if ((first < tf2.first) || (first == tf2.first && second < tf2.second) ||
        (first == tf2.first && second == tf2.second && third < tf2.third)) {
      return true;
    } else {
      return false;
    }
  }
  bool operator==(const TF& tf2) const
  {
    if (first == tf2.first && second == tf2.second && third == tf2.third) {
      return true;
    } else {
      return false;
    }
  }
  std::size_t operator()(const TF& tf) const
  {
    using std::hash;
    using std::size_t;

    return ((hash<int>()(tf.first) ^ (hash<int>()(tf.second) << 1)) >> 1) ^
           (hash<int>()(tf.third) << 1);
  }
};

class MeshWrapper
{
protected:
  struct GeometryData
  {
    std::vector<Eigen::Vector3d> vertex_coordinates_{};
    std::vector<Eigen::Vector3d> cell_centers_{};
  };

  std::unique_ptr<GeometryData> data_{ nullptr };

  /**
   *@brief check if  coordinates of cell centers have been calculated, if not
   *construct it
   */
  bool update_cell_centors();

  bool update_vertex_coordinates();

  /**
   * @brief precalulate some data for future use
   *
   * points store coordinates of vertices in Eigen form
   *
   */
  void setUp();

private:
  /*********** Functions begin **************/

  void readFromInp(std::ifstream& fin);

  void readFromOvm(std::ifstream& fin);

  /** insert a tetrahedral cell to mesh
      @param v  vertices of the tet
      @param faces  map from halfface(TF) to face, tell whether the face
                    it belong to had been added. faces shold be empty when
                    first use addCell.
  */
  void addCell(std::vector<OvmVeH>& v, std::map<TF, OvmFaH>& faces);

  /**
   *@brief given 4 vertices of a tet, push all 4 faces to the vector
   */
  void tetFaces(std::vector<Eigen::Matrix<long long, 3, 1>>& faces,
                long long v[4]);

  bool isSameHalfface(const std::vector<int>& f1, const std::vector<int>& f2);

public:
  /*********** Properties begin **************/

  std::unique_ptr<VMesh> ovm_mesh{ nullptr };

  std::string mesh_name = "";

  /*********** Properties end **************/

  // new ovm_mesh
  MeshWrapper()
    : data_{ std::make_unique<GeometryData>() }
    , ovm_mesh(new VMesh)
  {}

  // ovm_mesh is a shared_ptr
  ~MeshWrapper() {}

  /*********** IO related **************/

  void read_mesh(std::string filename);

  void save_as_OVM(std::string filename);

  // 分离文件的路径,文件名和后缀名
  // 0:path
  // 1:filename
  // 2:extension
  // 3:filename without extension
  std::vector<std::string> separate_file_name(std::string filename);

  /********************** geometry ****************************/

  std::string mesh_info() const;

  template<int VerticesNumber>
  void get_face_data(
    std::vector<Eigen::Vector3d>& points,
    std::vector<Eigen::Matrix<long long, VerticesNumber, 1>>& faces) const
  {
    size_t nv = ovm_mesh->n_vertices();
    points.reserve(nv);
    size_t nf = ovm_mesh->n_faces();
    faces.resize(nf);

    // get points data
    for (auto& p : this->data_->vertex_coordinates_) {
      points.push_back(p);
    }
    // get faces data
    for (auto fiter : ovm_mesh->faces()) {
      Eigen::Matrix<long long, VerticesNumber, 1> v;
      int k = 0;
      for (auto j = ovm_mesh->fv_iter(fiter); j.valid(); ++j) {
        v[k] = j->idx();
        k++;
      }
      faces[fiter.idx()] = v;
    }
  }

  /**
   *@brief return (points,faces)
   */
  std::tuple<RowMatrix<double>, RowMatrix<long long>> get_face_data() const;

  /**
   *@brief return (points,faces)
   */
  template<int VerticesNumber, typename Container>
  std::tuple<RowMatrix<double>, RowMatrix<long long>> get_face_data(
    const Container& face_ids) const
  {
    size_t nv = ovm_mesh->n_vertices();
    size_t nf = face_ids.size();
    RowMatrix<double> points(nv, VerticesNumber);
    RowMatrix<long long> faces(nf, VerticesNumber);

    // points
    for (size_t i = 0; i < nv; i++) {
      points.row(i) = data_->vertex_coordinates_[i];
    }

    // faces
    size_t i = 0;
    for (auto fi : face_ids) {
      int k = 0;
      OvmFaH fiter(fi);
      for (auto j = ovm_mesh->fv_iter(fiter); j.valid(); ++j) {
        faces(i, k) = j->idx();
        k++;
      }
      i++;
    }
    return std::tie(points, faces);
  }

  void get_boundary_faces(std::vector<int>& faceids_list) const;

  /**
  * @brief given cells' ids, return their faces's ids
  *
  * @param complement_set if true, return the cells' face ids that are not in
  cell_ids.
  */
  std::set<int> get_cell_face_ids(const std::set<int>& cell_ids,
                                  bool complement_set = false) const;

  /**
   *@brief return ids of boundary faces.
   */
  std::vector<int> get_boundary_face_ids();

  void get_separated_cells_data(
    std::vector<Eigen::Vector3d>& points,
    std::vector<Eigen::Matrix<long long, 3, 1>>& faces,
    std::vector<int>* ids_ptr = nullptr) const;

  void get_edge_data(std::vector<OvmEgH>& edge_ids,
                     std::vector<Eigen::Vector3d>& edge_points);
  /**
   *@reture average size(radius) of cells, for frame render
   */
  double cell_radius() const;

  double average_edge_length() const;

  /**
   * @brief return cell centors in matrix form
   */
  RowMatrix<double> request_cell_centers();

  /**
   *@brief return cells' volume
   */
  std::vector<double> request_cell_volumes();

  OvmFaH commonFace(OvmCeH ch1, OvmCeH ch2);

  const std::vector<Eigen::Vector3d>& get_vertex_coordinates()
  {
    if (data_->vertex_coordinates_.size() != ovm_mesh->n_vertices()) {
      update_vertex_coordinates();
    }
    return data_->vertex_coordinates_;
  };

  const std::vector<Eigen::Vector3d>& get_vertex_coordinates() const
  {
    return data_->vertex_coordinates_;
  }

  const std::vector<Eigen::Vector3d>& cell_centers()
  {
    if (data_->cell_centers_.size() != ovm_mesh->n_cells()) {
      update_cell_centors();
    }
    return data_->cell_centers_;
  }

  const std::vector<Eigen::Vector3d>& cell_centers() const
  {
    return data_->cell_centers_;
  }

  /******************* mesh related *************************/

  const VMesh* mesh_ptr() { return ovm_mesh.get(); }

  size_t n_vertices() const { return ovm_mesh->n_vertices(); }

  /**
   * @brief number of cells
   */
  size_t n_cells() const { return ovm_mesh->n_cells(); };

  /**
   * @brief find cells around an halfedge, ordered clockwise or
   anticlockwise
   * @warning make sure construct_matching_graph being run, and
   matching_graph *exits
   */
  std::set<size_t> get_adjacent_cells(OvmHaEgH hf_handle);
};