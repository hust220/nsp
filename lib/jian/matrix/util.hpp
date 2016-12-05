#pragma once

#include <iostream>
#include <fstream>
#include "traits.hpp"
#include "../utils/Serial.hpp"
#include "../utils/file.hpp"

BEGIN_JN


	template<typename NumType>
	using MatX = Eigen::Matrix<NumType, -1, -1>;
	using Mat = MatX<val_t>;
	using Matf = MatX<float>;
	using Matd = MatX<double>;
	using Mati = MatX<int>;
	using Matu = MatX<unsigned>;

	template<typename NumType>
	using VecX = Eigen::Matrix<NumType, -1, 1>;
	using Vec = VecX<val_t>;
	using Vecf = VecX<float>;
	using Vecd = VecX<double>;
	using Veci = VecX<int>;
	using Vecu = VecX<unsigned>;

	template<typename NumType>
	using RowVecX = Eigen::Matrix<NumType, 1, -1>;
	using RowVec = RowVecX<val_t>;
	using RowVecf = RowVecX<float>;
	using RowVecd = RowVecX<double>;
	using RowVeci = RowVecX<int>;
	using RowVecu = RowVecX<unsigned>;

	template<typename _NumType>
	inline std::ostream &operator <(std::ostream &stream, const MatX<_NumType> &mat) {
		stream < int(mat.rows()) < int(mat.cols());
		for (int i = 0; i < mat.rows(); i++) {
			for (int j = 0; j < mat.cols(); j++) {
				stream < mat(i, j);
			}
		}
		return stream;
	}

	template<typename _NumType>
	inline std::istream &operator >(std::istream &stream, MatX<_NumType> &mat) {
		int rows, cols;
		stream > rows > cols;
		mat.resize(rows, cols);
		for (int i = 0; i < mat.rows(); i++) {
			for (int j = 0; j < mat.cols(); j++) {
				stream > mat(i, j);
			}
		}
		return stream;
	}

	template<typename _Row>
	inline void mat_set_rows(Mat &mat, int beg, const _Row &row) {
		for (int i = 0; i < 3; i++) {
			mat(beg, i) = row[i];
		}
	}

	template<typename _First, typename _Second, typename... _Rest>
	inline void mat_set_rows(Mat &mat, int beg, const _First &first, const _Second &second, const _Rest &...rest) {
		mat_set_rows(mat, beg, first);
		mat_set_rows(mat, beg + 1, second, rest...);
	}

	template<typename _Vec, typename _Row>
	inline void vec_set(_Vec &vec, const _Row &row) {
		for (int i = 0; i < 3; i++) vec[i] = row[i];
	}

	template<typename _Vec, typename _Fn, typename _First, typename... _Rest>
	inline void vec_set(_Vec &vec, _Fn &&fn, const _First &first, const _Rest &...rest) {
		for (int i = 0; i < 3; i++) vec[i] = fn(first[i], (rest[i])...);
	}

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
	template<typename T = Eigen::MatrixXd, std::enable_if_t<is_stl_mat<T>::value, int> = 42>
	inline T make_mat(int row, int col) {
		using Vec = std::remove_reference_t<decltype(std::declval<T>().front())>;
		return T(row, Vec(col, 0));
	}

	template<typename T = Eigen::MatrixXd, std::enable_if_t<is_eigen_mat<T>::value, int> = 42>
	inline T make_mat(int row, int col) {
		return T::Zero(row, col);
	}

	// mat_from_file
	template<typename T = Eigen::MatrixXd>
	inline T mat_from_file(const Str &file) {
		std::ifstream ifile;
		FOPEN(ifile, file);
		int rows, cols;
		ifile >> rows >> cols;
		auto mat = make_mat<T>(rows, cols);
		for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) {
			if (ifile >> ref(mat, i, j)) continue; else throw "jian::mat_from_file failed!";
		}
		FCLOSE(ifile);
		return mat;
	}

	// rows_mats
	template<typename _Mat>
	inline int rows_mats(_Mat &&mat) {
		return mat.rows();
	}

	template<typename _First, typename _Second, typename... _Rest>
	inline int rows_mats(_First &&first, _Second &&second, _Rest && ...rest) {
		return rows_mats(first) + rows_mats(second, rest...);
	}

	// cols_mats
	inline int cols_mats() {
		return 0;
	}

	template<typename T, typename... L>
	inline int cols_mats(const T &mat, const L & ...mats) {
		return std::max(int(mat.cols()), cols_mats(mats...));
	}

	// hstack
	inline void hstack_helper(Mat &t, int n) {}

	template<typename _First, typename... _Rest>
	inline void hstack_helper(Mat &mat, int n, const _First &first, const _Rest & ...rest) {
		for (int i = 0; i < first.rows(); i++) {
			for (int j = 0; j < first.cols(); j++) mat(n, j) = first(i, j);
			n++;
		}
		hstack_helper(mat, n, rest...);
	}

	template<typename... _Mats>
	inline Mat hstack(const _Mats & ...mats) {
		Mat t(rows_mats(mats...), cols_mats(mats...));
		hstack_helper(t, 0, mats...);
		return t;
	}

END_JN

