#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <queue>
#include <set>
#include <tuple>
#include <vector>

template<typename T>
using SMatrix = Eigen::SparseMatrix<T, Eigen::RowMajor>;

template<typename T>
using Triple = Eigen::Triplet<T>;
using Tripled = Triple<double>;

using _Matrix_3 = double[3][3];

using Matrix_3 = Eigen::Matrix<double, 3, 3>;

using Triangle = Eigen::Matrix<long long, 3, 1>;

template<int n>
using FaceVertices = Eigen::Matrix<long long, n, 1>;

using V3d = Eigen::Vector3d;

using V4d = Eigen::Vector4d;

template<typename T>
using RowMatrix =
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
template<int columns>
using RowVectors =
  Eigen::Matrix<double, Eigen::Dynamic, columns, Eigen::RowMajor>;

using IndexType = long long;

const double PI = 3.141592653589793238463;
