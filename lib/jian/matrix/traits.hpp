#pragma once

#include <string>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Geometry>
#include <jian/utils/traits.hpp>
#include <vector>
#include <list>
#include <deque>
#include <array>

BEGIN_JN

template<typename T>
struct is_eigen_mat {
private:
	using F = std::decay_t<T>;
	template<typename _Scalar, int... Pars>
	static std::true_type check(Eigen::Matrix<_Scalar, Pars...>);
	static std::false_type check(...);
public:
	enum { value = std::is_same<decltype(check(std::declval<F>())), std::true_type>::value };
};

template<typename T>
struct is_stl_mat {
private:
	using F = std::decay_t<T>;
	template<template<typename...> class LS2, typename... Pars>
	static std::true_type check(std::vector<LS2<Pars...>>);
	template<template<typename...> class LS2, typename... Pars>
	static std::true_type check(std::deque<LS2<Pars...>>);
	template<template<typename...> class LS2, typename... Pars>
	static std::true_type check(std::list<LS2<Pars...>>);
	template<template<typename...> class LS2, typename... Pars, int N>
	static std::true_type check(std::array<LS2<Pars...>, N>);
	static std::false_type check(...);
public:
	enum { value = std::is_same<decltype(check(std::declval<F>())), std::true_type>::value };
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
	enum { value = std::is_same<decltype(check(std::declval<F>())), std::true_type>::value };
};

template<typename T>
struct mat_value_type {
private:
	template<template<typename...> class P, template<typename...> class U, typename F, typename... V,
		std::enable_if_t<is_stl_mat<T>::value, int> = 42>
		static F check(P<U<F, V...>>);
	template<typename U, int... F,
		std::enable_if_t<is_eigen_mat<T>::value, int> = 42>
		static U check(Eigen::Matrix<U, F...>);
public:
	using type = decltype(check(std::declval<T>()));
};

template<typename T>
using mat_value_type_t = typename mat_value_type<T>::type;

END_JN

