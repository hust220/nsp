#ifndef JIAN_ETL_MATRIX_TRAITS
#define JIAN_ETL_MATRIX_TRAITS

#include <cstring>
#include "../../Eigen/Dense"
#include "../../Eigen/SVD"
#include "../../Eigen/Geometry"
#include "../traits.h"
using namespace Eigen;

namespace jian {

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

} // namespace jian

#endif


