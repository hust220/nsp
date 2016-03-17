#ifndef JIAN_ETL_MATRIX_UTIL
#define JIAN_ETL_MATRIX_UTIL

#include "traits.hpp"

namespace jian {

// rows
template<typename MatType, std::enable_if_t<is_stl_mat<MatType>::value, int> = 42>
inline int rows(MatType &&mat) {
    return mat.size();
}

template<typename MatType, std::enable_if_t<is_eigen_mat<MatType>::value, int> = 42>
inline int rows(MatType &&mat) {
    return mat.rows();
}

template<typename MatType, std::enable_if_t<is_array_mat<MatType>::value, int> = 42>
inline int rows(MatType &&mat) {
    return std::distance(std::begin(mat), std::end(mat));
}

// ref
template<typename T, std::enable_if_t<is_stl_mat<T>::value, int> = 42>
inline template_first_parameter_t<template_first_parameter_t<T&&>> ref(T &&mat, int a, int b) {
    return mat[a][b];
}

template<typename T, std::enable_if_t<is_eigen_mat<T>::value, int> = 42>
inline template_first_parameter_t<T&&> ref(T &&mat, int a, int b) {
    return mat(a, b);
}

// make_mat
template<typename T = MatrixXd, std::enable_if_t<is_stl_mat<T>::value, int> = 42>
inline T make_mat(int row, int col) {
    using Vec = std::remove_reference_t<decltype(std::declval<T>().front())>;
    return T(row, Vec(col, 0));
}

template<typename T = MatrixXd, std::enable_if_t<is_eigen_mat<T>::value, int> = 42>
inline T make_mat(int row, int col) {
    return T::Zero(row, col);
}

// mat_from_file
template<typename T = MatrixXd>
inline T mat_from_file(const std::string &file) {
    std::ifstream ifile(file.c_str());
    int rows, cols;
    ifile >> rows >> cols;
    auto mat = make_mat<T>(rows, cols);
    for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) {
        if (ifile >> ref(mat, i, j)) continue; else throw "jian::mat_from_file failed!";
    }
    return mat;
}

// rows_mats
inline int rows_mats() {
    return 0;
}

template<typename T, typename... L>
inline int rows_mats(const T &mat, const L & ...mats) {
    return mat.rows() + rows_mats(mats...);
}

// cols_mats
inline int cols_mats() {
    return 0;
}

template<typename T, typename... L>
inline int cols_mats(const T &mat, const L & ...mats) {
    return std::max(mat.cols(), cols_mats(mats...));
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
    T t(rows_mats(mat, mats...), cols_mats(mat, mats...));
    hstack_helper(t, 0, mat, mats...);
    return t;
}

} // namespace jian

#endif


