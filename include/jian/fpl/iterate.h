#ifndef JIAN_FPL_ITERATE
#define JIAN_FPL_ITERATE

#include "../util/std.h"
#include "../etl/core.h"

namespace jian {
namespace fpl {

template<template<typename...> class ResultType = std::list, class Fn, class DataType, 
         template<typename...> class ListType, template<typename...> class... ListsType>
auto map(Fn &&f, const ListType<DataType> &list, const ListsType<DataType> & ...lists) 
        -> ResultType<decltype(f(list.front(), (lists.front())...))> {
    ResultType<decltype(f(list.front(), (lists.front())...))> result;
    for (int i = 0; i < std::distance(std::begin(list), std::end(list)); i++) {
        result.push_back(f(*(std::next(std::begin(list), i)), (*(std::next(std::begin(lists), i)))...));
    }
    return result;
}

//template<typename Fn, typename DataType, typename ListType, typename... ListsType, 
//         std::enable_if_t<etl::is_decay_same<std::pair<int, int>, ListType>::value and 
//                          etl::is_all_decay_same<std::pair<int, int>, ListsType...>::value, int> = 42>
//auto fold(Fn &&f, DataType &&data, ListType &&list, ListsType && ...lists) {
////         -> decltype(f(data, *std::begin(list), (*std::begin(lists))...)) {
//    using R = std::decay_t<decltype(f(data, list.first, (lists.first)...))>;
//    R sum = std::forward<DataType>(data);
//    for (int i = 0; i < list.second - list.first; i++) sum = f(sum, list.first + i, (lists.first + i)...);
//    return sum;
//}

template<typename Fn, typename DataType, typename ListType, typename... ListsType>
auto fold(Fn &&f, DataType &&data, ListType &&list, ListsType && ...lists) {
//         -> decltype(f(data, *std::begin(list), (*std::begin(lists))...)) {
    using R = std::decay_t<decltype(f(data, *std::begin(list), (*std::begin(lists))...))>;
    R sum = std::forward<DataType>(data);
    for (int i = 0; i < std::distance(std::begin(list), std::end(list)); i++) {
        sum = f(sum, *(std::next(std::begin(list), i)), (*(std::next(std::begin(lists), i)))...);
    }
    return sum;
}

template<typename Fn, typename ListType, typename... ListsType> 
void each(Fn &&f, ListType &&list, ListsType && ...lists) {
    int index = 0;
    for (int i = 0; i < std::distance(std::begin(list), std::end(list)); i++) {
        f(*(std::next(std::begin(list), i)), (*(std::next(std::begin(lists), i)))..., index);
        index++;
    }
}

template<typename Fn> 
void iterate2(const std::pair<int, int> &p1, const std::pair<int, int> &p2, Fn &&f) {
    for (int i = p1.first; i < p1.second; i++) {
        for (int j = p2.first; j < p2.second; j++) {
            f(i, j);
        }
    }
}

} // namespace fpl
} // namespace jian

#endif





