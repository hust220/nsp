#ifndef JIAN_ETL_CORE
#define JIAN_ETL_CORE

#include "../util/std.h"

namespace jian {
namespace etl {

template<typename T, typename U>
struct is_decay_same {
    enum {value = std::is_same<std::decay_t<T>, std::decay_t<U>>::value};
};

template<typename... T> struct is_all_decay_same;

template<typename T>
struct is_all_decay_same<T> {
    enum {value = true};
};

template<typename T, typename U, typename... F>
struct is_all_decay_same<T, U, F...> {
    enum {value = is_decay_same<T, U>::value and is_all_decay_same<U, F...>::value};
};

template<typename T, template<typename...> class U>
struct is_same_template {
private:
    template<template<typename...> class F, typename K, typename... Pars>
    K check(F<K, Pars...>);
public:
    enum {value = std::is_same<U<decltype(check(declval<T>()))>, std::decay_t<T>>::value};
};

template<typename T>
struct value_type {
private:
    template<template<typename...> class F, typename K, typename... Pars>
    K check(F<K, Pars...>);
public:
    using type = decltype(check(declval<T>()));
};

template<typename T>
using value_type_t = typename value_type<T>::type;

} // namespace etl
} // namespace jian

#endif





