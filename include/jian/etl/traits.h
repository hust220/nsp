#ifndef JIAN_ETL_TRAITS
#define JIAN_ETL_TRAITS

#include <type_traits>

namespace jian {

template<typename T, typename U>
struct is_decay_same {
    enum {value = std::is_same<std::decay_t<T>, std::decay_t<U>>::value};
};

template<typename... T>
struct is_all_decay_same;

template<typename T>
struct is_all_decay_same<T> {
    enum {value = true};
};

template<typename T, typename U, typename... F>
struct is_all_decay_same<T, U, F...> {
    enum {value = is_decay_same<T, U>::value && is_all_decay_same<U, F...>::value};
};

template<typename T, template<typename...> class U>
struct is_template_same {
private:
    template<template<typename...> class F, typename K, typename... Pars>
    static K check(F<K, Pars...>);
public:
    enum {value = std::is_same<U<decltype(check(std::declval<T>()))>, std::decay_t<T>>::value};
};

template<typename T>
struct template_first_parameter {
private:
    template<template<typename...> class F, typename K, typename... Pars>
    static K check(F<K, Pars...>);
    template<typename K, typename U> struct modify {using type = K;};
    template<typename K, typename U> struct modify<K, const U> {using type = const K;};
    template<typename K, typename U> struct modify<K, const U &> {using type = const K &;};
    template<typename K, typename U> struct modify<K, const U &&> {using type = const K &&;};
    template<typename K, typename U> struct modify<K, U&> {using type = K&;};
    template<typename K, typename U> struct modify<K, U&&> {using type = K&&;};
public:
    using type = typename modify<decltype(check(std::declval<T>())), T>::type;
};

template<typename T>
using template_first_parameter_t = typename template_first_parameter<T>::type;

} // namespace jian

#endif





