#ifndef JIAN_UTIL_MAT
#define JIAN_UTIL_MAT

#include <cstring>
#include "Point.h"
#include "Obj.h"
#include "MLib.h"
#include "../Eigen/Dense"
#include "../Eigen/SVD"
#include "../Eigen/Geometry"
using namespace Eigen;

namespace jian {
namespace mat {

// ref
template<typename MatType, typename DataType>
DataType &ref(MatType &mat, int a, int b);

template<typename MatType, typename DataType>
const DataType &ref(const MatType &mat, int a, int b);

template<typename DataType, int... Par>
auto & ref(Matrix<DataType, Par...> &mat, int a, int b) {
    return mat(a, b);
}

template<typename DataType, int... Par>
const auto & ref(const Matrix<DataType, Par...> &mat, int a, int b) {
    return mat(a, b);
}

template<typename MatType>
auto &ref(MatType &mat, int a, int b) {
    return mat[a][b];
}

template<typename MatType>
const auto &ref(const MatType &mat, int a, int b) {
    return mat[a][b];
}

// make_mat
MatrixXf make_mat(int row, int col) {
    return MatrixXf::Zero(row, col);
}

template<typename MatType>
MatType make_mat(int row, int col) {
    typedef typename std::remove_reference<decltype(MatType().front())>::type VecType;
    return MatType(row, VecType(col, 0));
}

template<>
MatrixXf make_mat<MatrixXf>(int row, int col) {
    return MatrixXf::Zero(row, col);
}

// mat_from_file
template<typename Mat = MatrixXf>
Mat mat_from_file(std::string file) {
    std::ifstream ifile(file.c_str());
    int rows, cols;
    ifile >> rows >> cols;
    Mat mat(rows, cols);
    for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) {
            if (ifile >> mat(i, j)) continue;
            else {
                std::cerr << "jian::mat_from_file failed!" << std::endl;
                exit(1);
            }
    }
    return mat;
}

// sum
template<typename F, typename T1>
auto sum(F &&f, T1 &&t1) ->decltype(f(t1)) {
    return f(std::forward<T1>(t1));
}

template<typename F, typename T1, typename... T2>
auto sum(F &&f, T1 &&t1, T2 && ...t2) ->decltype(f(t1)) {
    return f(std::forward<T1>(t1)) + sum(f, std::forward<T2>(t2)...);
}

// hstack
template<typename T1>
MatrixXf hstack(T1 &&mat1) {
    return std::forward<T1>(mat1);
}

template<typename T1, typename T2, typename... T3>
MatrixXf hstack(T1 &&mat1, T2 &&mat2, T3 && ...mat3) {
    if (mat1.rows() == 0) {
        return mat2;
    } else if (mat2.rows() == 0) {
        return mat1;
    }
    assert(mat1.cols() == mat2.cols());
    MatrixXf mat(mat1.rows() + mat2.rows(), mat1.cols());
    int i = 0;
    for (; i < mat1.rows(); i++) {
        for (int j = 0; j < mat1.cols(); j++) {
            mat(i, j) = mat1(i, j);
        }
    }
    for (int ii = 0; i < mat.rows(); i++, ii++) {
        for (int j = 0; j < mat2.cols(); j++) {
            mat(i, j) = mat2(ii, j);
        }
    }
    return hstack(mat, std::forward<T3>(mat3)...);
}

} // namespace mat
} // namespace jian

#endif


