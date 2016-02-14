#ifndef JIAN_UTIL_MAT
#define JIAN_UTIL_MAT

#include <cstring>
#include "Point.h"
#include "Obj.h"
#include "MLib.h"
#include "../Eigen/Dense"
#include "../Eigen/SVD"
#include "../Eigen/Geometry"
#include "../etl/control.h"
using namespace Eigen;

namespace jian {
namespace mat {

template<typename T>
struct is_eigen_mat {
private:
    using F = std::decay_t<T>;
    template<typename _Scalar, int... Pars> 
    static std::true_type check(Eigen::Matrix<_Scalar, Pars...>);
    static std::false_type check(...);
public:
    enum {value = std::is_same<decltype(check(declval<F>())), std::true_type>::value};
};

template<typename T>
struct is_stl_mat {
private:
    using F = std::decay_t<T>;
    template<template<typename...> class LS1, template<typename...> class LS2, typename... Pars> 
    static std::true_type check(LS1<LS2<Pars...>>);
    static std::false_type check(...);
public:
    enum {value = std::is_same<decltype(check(declval<F>())), std::true_type>::value};
};

template<typename T>
struct is_array_mat {
private:
    using F = std::decay_t<T>;
    template<typename U, std::size_t N> 
    static std::true_type check(U(*)[N]);
    template<typename U> 
    static std::true_type check(U **);
    static std::false_type check(...);
public:
    enum {value = std::is_same<decltype(check(declval<F>())), std::true_type>::value};
};

// rows
template<typename MatType, std::enable_if_t<is_eigen_mat<MatType>::value, int> = 42>
int rows(MatType &&mat) {
    return mat.rows();
}

template<typename MatType, std::enable_if_t<is_stl_mat<MatType>::value, int> = 42>
int rows(MatType &&mat) {
    return mat.size();
}

template<typename MatType, std::enable_if_t<is_array_mat<MatType>::value, int> = 42>
int rows(MatType &&mat) {
    return std::distance(std::begin(mat), std::end(mat));
}

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
MatrixXd make_mat(int row, int col) {
    return MatrixXd::Zero(row, col);
}

template<typename MatType>
MatType make_mat(int row, int col) {
    typedef typename std::remove_reference<decltype(MatType().front())>::type VecType;
    return MatType(row, VecType(col, 0));
}

template<>
MatrixXd make_mat<MatrixXd>(int row, int col) {
    return MatrixXd::Zero(row, col);
}

// mat_from_file
template<typename Mat = MatrixXd>
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

// rows_mats
inline int rows_mats() {
    return 0;
}

template<typename T, typename... L>
inline int rows_mats(const T &mat, const L & ...mats) {
    return mat.rows() + rows_mats(mats...);
}

// hstack
template<typename T>
inline void hstack_helper(T &t, int n) {}

template<typename T, typename... L>
inline void hstack_helper(T &t, int n, const T &mat, const L & ...mats) {
    for (int i = 0; i < mat.rows(); i++) {
        for (int j = 0; j < mat.cols(); j++) t(n, j) = mat(i, j);
        n++;
    }
    hstack_helper(t, n, mats...);
}

template<typename T, typename... L>
inline T hstack(const T &mat, const L & ...mats) {
    T t(rows_mats(mat, mats...), mat.cols());
    hstack_helper(t, 0, mat, mats...);
    return t;
}

} // namespace mat
} // namespace jian

#endif


