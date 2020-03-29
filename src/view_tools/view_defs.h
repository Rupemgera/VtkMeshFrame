/**
 * @file view_defs.h
 * @author Liao Yizhou (3130001008@zju.edu.cn)
 * @brief alias, class for data exchange
 * @version 0.1
 * @date 2020-01-03
 *
 * @copyright Copyright (c) 2020
 *
 */

#pragma once
#include <Eigen/Dense>

/***************** defines ********************/

using Point3d = Eigen::Vector3d;
template <int n>
// VertexList<n> is a n-list storing the vertices of a facet
using VertexList = Eigen::Matrix<long long, n, 1>;
using IndexType = long long;
template <typename T>
using RowMatrix =
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

/***************** classes ********************/

class Color {
public:
  Color(double *data) {
    for (size_t i = 0; i < 3; i++) {
      _data[i] = data[i];
    }
  }
  Color(double r, double g, double b) {
    _data[0] = r;
    _data[1] = g;
    _data[2] = b;
  }
  Color(std::string description) {
    if (description == "red") {
      _data[0] = 1.0;
      _data[1] = 0.0;
      _data[2] = 0.0;
    } else if (description == "green") {
      _data[0] = 0.0;
      _data[1] = 1.0;
      _data[2] = 0.0;
    } else if (description == "blue") {
      _data[0] = 0.0;
      _data[1] = 0.0;
      _data[2] = 1.0;
    } else {
      _data[0] = 1.0;
      _data[1] = 1.0;
      _data[2] = 1.0;
    }
  }
  const double *data() const { return _data; }
  double *data() { return _data; }
  const double operator[](size_t n) const { return _data[n % 3]; }

private:
  double _data[3];
};

class VtkWrapperConfiguration {
public:
  struct MeshConfig {
    int edge_color[3];
    int face_color[3];
  };
};