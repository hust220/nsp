#pragma once

#include <list>
#include <vector>
#include <algorithm>
#include <deque>
#include "traits.hpp"

BEGIN_JN

template<typename LS>
inline int len(LS &&ls) {
    return std::distance(std::begin(ls), std::end(ls));
}

template<typename LS, std::enable_if_t<is_template_same<LS, std::list>::value, int> = 42>
inline template_first_parameter_t<LS&&> ref(LS &&ls, int n) {
    if (n < 0) return *std::next(ls.end(), n);
    else return *std::next(ls.begin(), n);
}

template<typename LS, std::enable_if_t<!is_template_same<LS, std::list>::value, int> = 42>
inline template_first_parameter_t<LS&&> ref(LS &&ls, int n) {
    if (n < 0) return ls[len(ls) + n];
    else return ls[n];
}

template<typename Fn, typename ListType, typename... ListsType>
bool exists(Fn &&f, const ListType &list, const ListsType & ...lists) {
    for (int i = 0; i < list.size(); i++) {
        if (f(*(std::next(std::begin(list), i)), (*(std::next(std::begin(lists), i)))...)) return true;
    }
    return false;
}

template<typename Fn, typename ListType>
auto find(Fn &&f, const ListType &list) {
    for (int i = 0; i < list.size(); i++) {
        const auto &val = *std::next(std::begin(list), i);
        if (f(val)) return val;
    }
    throw "not found!";
}

template<typename Fn> 
void iterate2(std::pair<int, int> p1, std::pair<int, int> p2, Fn &&f) {
    for (int i = p1.first; i < p1.second; i++) {
        for (int j = p2.first; j < p2.second; j++) {
            f(i, j);
        }
    }
}

template<template<typename...> class ResultType = std::list, class Fn, class DataType, template<typename...> class ListType = std::vector>
ResultType<DataType> filter(Fn &&f, const ListType<DataType> &list) {
    ResultType<DataType> result;
    for (auto && data : list) if (f(data)) result.push_back(data);
    return result;
}

template<template<typename...> class ResultType = std::list, class Fn, class DataType, 
         template<typename...> class ListType, template<typename...> class... ListsType>
auto map(Fn &&f, const ListType<DataType> &list, const ListsType<DataType> & ...lists) -> ResultType<decltype(f(list.front(), (lists.front())...))>{
    ResultType<decltype(f(list.front(), (lists.front())...))> result;
    for (int i = 0; i < std::distance(std::begin(list), std::end(list)); i++) {
        result.push_back(f(*(std::next(std::begin(list), i)), (*(std::next(std::begin(lists), i)))...));
    }
    return result;
}

template<typename Fn, typename ListType, typename... ListsType> 
void each(Fn &&f, ListType &&list, ListsType && ...lists) {
    int index = 0;
    for (int i = 0; i < std::distance(std::begin(list), std::end(list)); i++) {
        f(*(std::next(std::begin(list), i)), (*(std::next(std::begin(lists), i)))..., index);
        index++;
    }
}

template<typename Fn, typename DataType>
auto fold(Fn &&f, DataType &&data, std::pair<int, int> p) -> decltype(f(data, 0)) {
    decltype(f(data, 0)) sum = data;
    for (int i = p.first; i < p.second; i++) sum = f(sum, i);
    return sum;
}

template<typename Fn, typename DataType, typename ListType, typename... ListsType>
auto fold(Fn &&f, DataType &&data, ListType &&list, ListsType && ...lists) -> decltype(f(data, *std::begin(list), (*std::begin(lists))...)) {
    decltype(f(data, *std::begin(list), (*std::begin(lists))...)) sum = std::forward<DataType>(data);
    for (int i = 0; i < std::distance(std::begin(list), std::end(list)); i++) {
        sum = f(sum, *(std::next(std::begin(list), i)), (*(std::next(std::begin(lists), i)))...);
    }
    return sum;
}

template<typename T = std::deque<int>>
T range(double beg, double end, double step = 1) {
    using val_t = template_first_parameter_t<T>;
    T ls; for (val_t i = beg; i < end; i += step) ls.push_back(i);
    return ls;
}

template<typename List> void append(List &list) {}

template<typename List, typename Data1, typename... Data2> 
void append(List &list, Data1 &&data, Data2 && ...datum) {
    list.push_back(std::forward<Data1>(data));
    append(list, std::forward<Data2>(datum)...);
}

template<typename List> 
List slice(const List &list, int m, int n) {
    List new_list;
    std::copy(std::next(std::begin(list), m), std::next(std::begin(list), n), std::back_inserter(new_list));
    return new_list;
}

END_JN

