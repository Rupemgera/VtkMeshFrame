#include "mesh_wrapper.h"
#include <assert.h>
#include <regex>
#include <sstream>

void
MeshWrapper::readFromInp(std::ifstream& fin)
{
  /* 读入.inp文件 */

  std::string line_str;

  /* 正则表达式  */

  //匹配 *NODE信息
  std::regex node_regex("^\\*node", std::regex::icase);

  //匹配 *ELEMENT信息
  std::regex element_regex("^\\*element", std::regex::icase);

  //匹配 *
  std::regex comments_regex("\\*{2}|^$");

  //输入格式

  std::regex node_format(
    "^([\\d\\.Ee\\s]+),([\\d\\.Ee\\s]+),([\\d\\.Ee\\s]+),([\\d\\.Ee\\s]+)$");
  std::regex element_format("^([\\d\\.Ee\\s]+),([\\d\\.Ee\\s]+),([\\d\\.Ee\\s]+"
                            "),([\\d\\.Ee\\s]+),([\\d\\.Ee\\s]+)$");

  /*  读入点信息  */

  std::map<TF, OvmFaH> faces;
  faces.clear();

  std::string tmp1, tmp2;
  bool node_section = false, element_section = false;

  while (std::getline(fin, line_str)) {
    if (line_str == "")
      continue;
    if (line_str[0] == '*') {
      if (std::regex_search(line_str, node_regex)) {
        node_section = true;
        element_section = false;
      } else if (std::regex_search(line_str, element_regex)) {
        node_section = false;
        element_section = true;
      }
    } else {
      if (node_section) {
        //匹配到*NODE，开始从tmp1读入点信息
        int id;
        double x, y, z;
        // std::cout << tmp1.c_str() <<std::endl;
#ifdef __linux
        sscanf(line_str.c_str(), "%d, %lf , %lf , %lf", &id, &x, &y, &z);
#endif // __linux
#ifdef _WIN64
        sscanf_s(line_str.c_str(), "%d, %lf , %lf , %lf", &id, &x, &y, &z);
#endif // _WIN64 \
	//插入网格点
        ovm_mesh->add_vertex(MeshPoint(x, y, z));
      } else if (element_section) {
        int id;
        int v[4];
#ifdef __linux
        sscanf(line_str.c_str(),
               "%d, %d ,%d , %d , %d",
               &id,
               v,
               v + 1,
               v + 2,
               v + 3);
#endif // __linux
#ifdef _WIN64
        sscanf_s(line_str.c_str(),
                 "%d, %d ,%d , %d , %d",
                 &id,
                 v,
                 v + 1,
                 v + 2,
                 v + 3);
#endif // _WIN64 \
	//int转换为handle
        std::vector<OvmVeH> vs;
        vs.clear();
        // inp index start from 1 while OVM from 0
        vs.push_back((OvmVeH)(v[0] - 1));
        vs.push_back((OvmVeH)(v[1] - 1));
        vs.push_back((OvmVeH)(v[2] - 1));
        vs.push_back((OvmVeH)(v[3] - 1));
        // add cell
        addCell(vs, faces);
      }
    }
  }
}

void
MeshWrapper::readFromOvm(std::ifstream& fin)
{
  OpenVolumeMesh::IO::FileManager manager;
  // only print error infomation
  manager.setVerbosityLevel(1);
  manager.readStream(fin, *ovm_mesh);
}

void
MeshWrapper::save_as_OVM(std::string filename)
{
  OpenVolumeMesh::IO::FileManager manager;
  // OpenVolumeMesh::GeometricPolyhedralMeshV3d *mesh;
  manager.writeFile(filename, *ovm_mesh);
}

std::vector<std::string>
MeshWrapper::separate_file_name(std::string filename)
{
  std::vector<std::string> ret;
  ret.clear();

  size_t slash_position = filename.find_last_of('/');
  if (slash_position == filename.npos) {
    slash_position = 0;
  }
  size_t dot_position = filename.find_last_of('.');

  /********** path *************/
  ret.push_back(filename.substr(0, slash_position));
  /********** filename *************/
  ret.push_back(filename.substr(slash_position + 1));
  /********** extension *************/
  ret.push_back(filename.substr(dot_position + 1));
  /********** filename without extension *************/
  ret.push_back(filename.substr(0, dot_position));

  return ret;
}

std::string
MeshWrapper::mesh_info() const
{
  // std::string info = "";
  std::stringstream ss;

  ss << "vertices : ";
  ss << ovm_mesh->n_vertices();
  ss << '\n';
  ss << "edges : ";
  ss << ovm_mesh->n_edges();
  ss << '\n';
  ss << "faces : ";
  ss << ovm_mesh->n_faces();
  ss << '\n';
  ss << "cells : ";
  ss << ovm_mesh->n_cells();
  ss << '\n';
  // ss >> info;

  return ss.str();
}

std::tuple<RowMatrix<double>, RowMatrix<long long>>
MeshWrapper::get_face_data() const
{
  size_t nv = ovm_mesh->n_vertices();
  size_t nf = ovm_mesh->n_faces();
  RowMatrix<double> points(nv, 3);
  RowMatrix<long long> faces(nf, 3);
  // points
  for (auto viter : ovm_mesh->vertices()) {
    auto p = ovm_mesh->vertex(viter);
    size_t i = viter.idx();
    for (size_t j = 0; j < 3; j++) {
      points(i, j) = p[j];
    }
  }

  // faces
  for (auto fiter : ovm_mesh->faces()) {
    int k = 0;
    size_t i = fiter.idx();
    for (auto j = ovm_mesh->fv_iter(fiter); j.valid(); ++j) {
      faces(i, k) = j->idx();
      k++;
    }
  }

  return std::tie(points, faces);
}

void
MeshWrapper::get_boundary_faces(std::vector<int>& faceids_list) const
{
  faceids_list.clear();
  for (auto fiter : ovm_mesh->faces()) {
    if (ovm_mesh->is_boundary(fiter)) {
      faceids_list.push_back(fiter.idx());
    }
  }
}

std::set<int>
MeshWrapper::get_cell_face_ids(const std::set<int>& cell_ids,
                               bool complement_set) const
{
  std::set<int> ids;

  if (complement_set) {
    size_t n = ovm_mesh->n_cells();
    for (int i = 0; i < n; i++) {
      if (cell_ids.find(i) == cell_ids.end()) {
        auto cfiter = ovm_mesh->cf_iter(static_cast<OvmCeH>(i));
        for (; cfiter.valid(); ++cfiter) {
          ids.insert(cfiter->idx());
        }
      }
    }
  } else {
    for (int c : cell_ids) {
      auto cfiter = ovm_mesh->cf_iter(static_cast<OvmCeH>(c));
      for (; cfiter.valid(); ++cfiter) {
        ids.insert(cfiter->idx());
      }
    }
  }

  return ids;
}

std::vector<int>
MeshWrapper::get_boundary_face_ids()
{
  std::vector<int> faceids_list;
  for (auto fiter : ovm_mesh->faces()) {
    if (ovm_mesh->is_boundary(fiter)) {
      faceids_list.push_back(fiter.idx());
    }
  }
  return faceids_list;
}

void
MeshWrapper::get_separated_cells_data(
  std::vector<Eigen::Vector3d>& points,
  std::vector<Eigen::Matrix<long long, 3, 1>>& faces,
  std::vector<int>* ids_ptr) const
{
  // function to get geometric info of shrinked cells
  auto get_shrinked_data = [this, &points, &faces](OvmCeH& citer) {
    double rate = 0.2;
    int k = (int)points.size();
    std::vector<OvmVec3d> vertices;
    OvmVec3d tet_center(0.0, 0.0, 0.0);

    // calculate center point of tet
    for (auto vciter = ovm_mesh->cv_iter(citer); vciter.valid(); ++vciter) {
      vertices.push_back(ovm_mesh->vertex(*vciter));
      tet_center += ovm_mesh->vertex(*vciter);
    }
    tet_center /= vertices.size();
    // every vertices move a step to center
    for (auto& v : vertices) {
      v = v + (tet_center - v) * rate;
      points.push_back(Eigen::Vector3d(v.data()));
    }

    long long v[4] = { k, k + 1, k + 2, k + 3 };
    // 将四个顶点构成的tet的四个面push进faces
    for (int i = 0; i < 4; ++i) {
      Eigen::Matrix<long long, 3, 1> face_vertices;
      int q = 0;
      for (auto j = 0; j < 4; ++j) {
        if (j != i) {
          face_vertices[q] = v[j];
          q++;
        }
      }
      faces.push_back(face_vertices);
    }
  };

  size_t n;
  if (ids_ptr == nullptr) {
    n = ovm_mesh->n_cells();
    faces.clear();
    points.clear();
    faces.reserve(n * 4);
    points.reserve(n * 4);
    for (auto citer : ovm_mesh->cells()) {
      get_shrinked_data(citer);
    }
  } else {
    n = ids_ptr->size();
    faces.clear();
    points.clear();
    faces.reserve(n * 4);
    points.reserve(n * 4);
    for (size_t i = 0; i < n; i++) {
      OvmCeH citer((int)i);
      get_shrinked_data(citer);
    }
  }

  assert(faces.size() == n * 4);
  assert(points.size() == n * 4);
}

void
MeshWrapper::get_edge_data(std::vector<OvmEgH>& edge_ids,
                           std::vector<Eigen::Vector3d>& edge_points)
{
  edge_points.reserve(edge_ids.size() * 2);
  for (auto eh : edge_ids) {
    OvmVeH u = ovm_mesh->edge(eh).to_vertex();
    OvmVeH v = ovm_mesh->edge(eh).from_vertex();
    edge_points.push_back(Eigen::Vector3d(ovm_mesh->vertex(u).data()));
    edge_points.push_back(Eigen::Vector3d(ovm_mesh->vertex(v).data()));
  }
}

double
MeshWrapper::cell_radius() const
{
  double maxr = 0, minr = 1e20;
  double sum = 0.0;
  // size is defined as  body center to vertex distance
  for (auto citer = ovm_mesh->cells_begin(); citer != ovm_mesh->cells_end();
       ++citer) {
    MeshPoint c(0, 0, 0);
    MeshPoint p[4];
    int k = 0;
    for (auto viter = ovm_mesh->cv_iter(*citer); viter.valid(); ++viter) {
      p[k] = ovm_mesh->vertex(*viter);
      c += p[k];
      k++;
    }
    c = c / 4;
    for (size_t i = 0; i < 4; i++) {
      sum += (p[i] - c).norm();
    }
  }
  return sum / ovm_mesh->n_cells();
}

double
MeshWrapper::average_edge_length() const
{
  double sum = 0.0;
  for (auto eit = ovm_mesh->halfedges_begin(); eit != ovm_mesh->halfedges_end();
       ++eit) {
    auto u = ovm_mesh->vertex(ovm_mesh->to_vertex_handle(*eit));
    auto v = ovm_mesh->vertex(ovm_mesh->from_vertex_handle(*eit));
    sum += (u - v).norm();
  }
  return sum / ovm_mesh->n_halfedges();
}

OvmFaH
MeshWrapper::commonFace(OvmCeH ch1, OvmCeH ch2)
{
  OvmFaH com_fh = ovm_mesh->InvalidFaceHandle;
  auto hfhs_vec1 = ovm_mesh->cell(ch1).halffaces(),
       hfhs_vec2 = ovm_mesh->cell(ch2).halffaces();

  std::set<OvmFaH> fhs_set1, fhs_set2;
  std::for_each(hfhs_vec1.begin(), hfhs_vec1.end(), [&](OvmHaFaH hfh) {
    fhs_set1.insert(ovm_mesh->face_handle(hfh));
  });
  std::for_each(hfhs_vec2.begin(), hfhs_vec2.end(), [&](OvmHaFaH hfh) {
    fhs_set2.insert(ovm_mesh->face_handle(hfh));
  });

  std::vector<OvmFaH> com_fhs;
  std::set_intersection(fhs_set1.begin(),
                        fhs_set1.end(),
                        fhs_set2.begin(),
                        fhs_set2.end(),
                        std::back_inserter(com_fhs));
  if (com_fhs.size() == 1)
    com_fh = com_fhs.front();
  return com_fh;
}

std::set<size_t>
MeshWrapper::get_adjacent_cells(OvmHaEgH hf_handle)
{
  std::set<size_t> res;
  auto citer = ovm_mesh->hec_iter(hf_handle);
  for (; citer.valid(); ++citer) {
    res.insert(citer->idx());
  }
  return res;
}

bool
MeshWrapper::update_cell_centors()
{
  // if property is not calculated
  data_->cell_centers_.resize(ovm_mesh->n_cells());
  for (auto citer : ovm_mesh->cells()) {
    MeshPoint c(0, 0, 0);
    for (auto viter = ovm_mesh->cv_iter(citer); viter.valid(); ++viter) {
      c += ovm_mesh->vertex(*viter);
    }
    data_->cell_centers_[citer.idx()] = Eigen::Vector3d((c / 4).data());
  }

  return true;
}

bool
MeshWrapper::update_vertex_coordinates()
{
  size_t n = ovm_mesh->n_vertices();
  data_->vertex_coordinates_.clear();
  data_->vertex_coordinates_.reserve(n);
  for (auto& v : ovm_mesh->vertices()) {
    data_->vertex_coordinates_.push_back(V3d(ovm_mesh->vertex(v).data()));
  }
  return true;
}

RowMatrix<double>
MeshWrapper::request_cell_centers()
{
  update_cell_centors();

  size_t n = data_->cell_centers_.size();
  RowMatrix<double> res(n, 3);

  for (size_t i = 0; i < n; ++i) {
    res.row(i) = data_->cell_centers_[i];
  }

  return res;
}

std::vector<double>
MeshWrapper::request_cell_volumes()
{
  size_t n = ovm_mesh->n_cells();
  std::vector<double> volumes(n);
  for (auto cit : ovm_mesh->cells()) {
    Eigen::Matrix<double, 4, 4> det;
    /*
                      | 1     1     1     1    |

       v=1/6 * det    | x1    x2   x3    x4 |

                      | y1   y2   y3   y4  |

                      | z1   z2    z3   z4  |
     */
    int k = 0;
    for (auto cvit = ovm_mesh->cv_iter(cit); cvit.valid(); ++cvit) {
      auto p = ovm_mesh->vertex(*cvit);
      det.col(k) << 1, p[0], p[1], p[2];
      k++;
    }
    volumes[cit.idx()] = abs(det.determinant()) / 6;
  }
  return volumes;
}

bool
MeshWrapper::isSameHalfface(const std::vector<int>& f1,
                            const std::vector<int>& f2)
{
  if ((f1[0] == f2[0] && f1[1] == f2[1]) ||
      (f1[1] == f2[0] && f1[2] == f2[1]) || (f1[2] == f2[0] && f1[0] == f2[1]))
    return true;
  return false;
}

void
MeshWrapper::setUp()
{
  // store vertices' coordinates in advance
  update_vertex_coordinates();
  // get cell centers
  update_cell_centors();
}

/**
 *添加面时顶点顺序为右手法则法向朝外
 */
void
MeshWrapper::addCell(std::vector<OvmVeH>& v, std::map<TF, OvmFaH>& faces)
{
  OvmFaH f0, f1, f2, f3;
  std::vector<OvmHaFaH> hf;
  std::vector<OvmVeH> fv;

  // f0 : 0->1->2
  TF tf0(v[0].idx(), v[1].idx(), v[2].idx());
  fv.clear();
  fv.push_back(v[0]);
  fv.push_back(v[1]);
  fv.push_back(v[2]);
  if (faces.find(tf0) == faces.end()) {
    f0 = ovm_mesh->add_face(fv);
    faces[tf0] = f0;
    hf.push_back(ovm_mesh->halfface_handle(f0, 0));
  } else {
    f0 = faces.find(tf0)->second;
    hf.push_back(ovm_mesh->halfface_handle(f0, 1));
  }

  // f1 : 0->2->3
  TF tf1(v[0].idx(), v[2].idx(), v[3].idx());
  fv.clear();
  fv.push_back(v[0]);
  fv.push_back(v[2]);
  fv.push_back(v[3]);
  if (faces.find(tf1) == faces.end()) {
    f1 = ovm_mesh->add_face(fv);
    faces[tf1] = f1;
    hf.push_back(ovm_mesh->halfface_handle(f1, 0));
  } else {
    f1 = faces.find(tf1)->second;
    hf.push_back(ovm_mesh->halfface_handle(f1, 1));
  }

  // f2 : 0->3->1
  TF tf2(v[0].idx(), v[3].idx(), v[1].idx());
  fv.clear();
  fv.push_back(v[0]);
  fv.push_back(v[3]);
  fv.push_back(v[1]);
  if (faces.find(tf2) == faces.end()) {
    f2 = ovm_mesh->add_face(fv);
    faces[tf2] = f2;
    hf.push_back(ovm_mesh->halfface_handle(f2, 0));
  } else {
    f2 = faces.find(tf2)->second;
    hf.push_back(ovm_mesh->halfface_handle(f2, 1));
  }

  // f3 : 1->3->2
  TF tf3(v[1].idx(), v[3].idx(), v[2].idx());
  fv.clear();
  fv.push_back(v[1]);
  fv.push_back(v[3]);
  fv.push_back(v[2]);
  if (faces.find(tf3) == faces.end()) {
    f3 = ovm_mesh->add_face(fv);
    faces[tf3] = f3;
    hf.push_back(ovm_mesh->halfface_handle(f3, 0));
  } else {
    f3 = faces.find(tf3)->second;
    hf.push_back(ovm_mesh->halfface_handle(f3, 1));
  }

  // add cell
  ovm_mesh->add_cell(hf);
}

void
MeshWrapper::tetFaces(std::vector<Eigen::Matrix<long long, 3, 1>>& faces,
                      long long v[4])
{
  for (int i = 0; i < 4; ++i) {
    Eigen::Matrix<long long, 3, 1> face_vertices;
    int k = 0;
    for (auto j = 0; j < 4; ++j) {
      if (j != i) {
        face_vertices[k] = v[j];
        k++;
      }
    }
    faces.push_back(face_vertices);
  }
}

void
MeshWrapper::read_mesh(std::string filename)
{
  std::ifstream fin(filename);
  if (!fin) {
    std::cerr << "file doesn't exist" << std::endl;
    return;
  }

  // 若之前已经载入过mesh了，则将之前的注销重新new一个空mesh
  auto names = separate_file_name(filename);
  mesh_name = names[3];
  // for (size_t i = 0; i < names.size(); i++)
  //   std::cerr << "name[" << i << "] : " << names[i] << std::endl;

  if (names[2] == "inp") {
    readFromInp(fin);
    std::string ovm_file_name = names[3] + ".ovm";
    save_as_OVM(ovm_file_name);
  } else if (names[2] == "ovm") {
    readFromOvm(fin);
  }
  std::cerr << "load succeed" << std::endl;

  // set up
  // setUp();
  update_cell_centors();
  update_vertex_coordinates();
}