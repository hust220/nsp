#ifndef JIAN_FPL_ITERATE
#define JIAN_FPL_ITERATE

#include "../util/std.h"

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

template<typename Fn, typename DataType>
auto fold(Fn &&f, DataType &&data, std::pair<int, int> p) -> decltype(f(data, 0)) {
    decltype(f(data, 0)) sum = data;
    for (int i = p.first; i < p.second; i++) sum = f(sum, i);
    return sum;
}

template<typename Fn, typename DataType, typename ListType, typename... ListsType>
auto fold(Fn &&f, DataType &&data, ListType &&list, ListsType && ...lists) 
         -> decltype(f(data, *std::begin(list), (*std::begin(lists))...)) {
    decltype(f(data, *std::begin(list), (*std::begin(lists))...)) sum = std::forward<DataType>(data);
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





