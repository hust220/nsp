#ifndef JIAN_MLIB_H
#define JIAN_MLIB_H

#include "std.h"
#include "mat.h"

namespace jian {

template<typename T, size_t N> char (&_ArraySizeHelper(T (&array)[N]))[N];

template<typename Array> std::size_t count(Array array) {
    return sizeof _ArraySizeHelper(array);
}

template<typename Fn>
auto let(Fn &&f) {
    return f();
}

template<template<typename...> class ListType, typename DataType>
DataType &ref(ListType<DataType> &list, int n);

template<template<typename...> class ListType, typename DataType>
const DataType &ref(const ListType<DataType> &list, int n);

template<typename DataType>
DataType &ref(std::list<DataType> &list, int n) {
    return *std::next(list.begin(), n);
}

template<typename DataType>
const DataType &ref(const std::list<DataType> &list, int n) {
    return *std::next(list.begin(), n);
}

template<template<typename...> class ListType, typename DataType>
DataType &ref(ListType<DataType> &list, int n) {
    return list[n];
}

template<template<typename...> class ListType, typename DataType>
const DataType &ref(const ListType<DataType> &list, int n) {
    return list[n];
}

template<typename DataType> DataType square(const DataType &data) {
    return data * data;
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

template<typename Fn> void iterate2(std::pair<int, int> p1, std::pair<int, int> p2, Fn &&f) {
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

template<typename Fn, typename ListType, typename... ListsType> void each(Fn &&f, ListType &&list, ListsType && ...lists) {
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

//template<typename Fn, typename DataType, typename ListType>
//auto fold(Fn &&f, DataType &&data, ListType &&list) -> decltype(f(data, *std::begin(list))) {
//    decltype(f(data, *std::begin(list))) sum = std::forward<DataType>(data);
//    for (auto && value : list) sum = f(sum, value);
//    return sum;
//}
//
template<typename Fn, typename DataType, typename ListType, typename... ListsType>
auto fold(Fn &&f, DataType &&data, ListType &&list, ListsType && ...lists) -> decltype(f(data, *std::begin(list), (*std::begin(lists))...)) {
    decltype(f(data, *std::begin(list), (*std::begin(lists))...)) sum = std::forward<DataType>(data);
    for (int i = 0; i < std::distance(std::begin(list), std::end(list)); i++) {
        sum = f(sum, *(std::next(std::begin(list), i)), (*(std::next(std::begin(lists), i)))...);
    }
    return sum;
}

template<typename List = std::vector<int>, typename DataType, typename StepType = int> List range(DataType &&beg, int len, StepType step = 1) {
    List list;
    list.reserve(len);
    DataType value = beg;
    for (int i = 0; i < len; i++) {
        list.push_back(value);
        value = value + step;
    }
    return list;
}

template<typename List> void append(List &list) {}

template<typename List, typename Data1, typename... Data2> void append(List &list, Data1 &&data, Data2 && ...datum) {
    list.push_back(std::forward<Data1>(data));
    append(list, std::forward<Data2>(datum)...);
}

template<typename List> List slice(const List &list, int m, int n) {
    List new_list;
    std::copy(std::next(std::begin(list), m), std::next(std::begin(list), n), std::back_inserter(new_list));
    return new_list;
}

inline void tokenize(const std::string &str, std::vector<string> &tokens, const std::string &delimiters) {
    tokens.clear();
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while (string::npos != pos || string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

inline void tokenize(const string &str, vector<string> &tokens, const string &delimiters, const string &temp) {
    tokens.clear();
    std::vector<std::pair<string::size_type, string::size_type>> vec;
    string::size_type first_i, first_j, second_i, second_j;
    int expected = 0;
    for (string::size_type i = 0; i < str.size(); i++) {
        int flag = 0;
        string::size_type j;
        for (j = 0; j < temp.size(); j++) {
            if (str[i] == temp[j]) {
                if (j % 2 == 0 && expected == 0) {
                    flag = 1;
                    break;
                } else if (j % 2 == 1 && expected == 1) {
                    flag = 2;
                    break;
                }
            }
        }
        if (flag == 1) {
            first_i = i;
            first_j = j;
            expected = 1;
        } else if (flag == 2 && j - first_j == 1) {
            second_i = i;
            second_j = j;
            expected = 0;
            vec.push_back(std::make_pair(first_i, second_i));
        }
    }
//    for (int i = 0; i < temp.size(); i += 2) {
//        string::size_type pos1 = str.find_first_of(string() + temp[i], 0);
//        string::size_type pos2 = str.find_first_of(string() + temp[i + 1], pos1 + 1);
//        while (string::npos != pos1 && string::npos != pos2) {
//            vec.push_back(make_pair(pos1, pos2));
//            pos1 = str.find_first_of(string() + temp[i], pos2 + 1);
//            pos2 = str.find_first_of(string() + temp[i + 1], pos1 + 1);
//        }
//    }
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while (any_of(vec.begin(), vec.end(), [&pos](const pair<string::size_type, string::size_type> &p){
        return pos != string::npos && p.first < pos && pos < p.second;
    })) {
        pos = str.find_first_of(delimiters, pos + 1);
    }
    while (string::npos != pos || string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
        while (any_of(vec.begin(), vec.end(), [&pos](const pair<string::size_type, string::size_type> &p){
            return pos != string::npos && p.first < pos && pos < p.second;
        })) {
            pos = str.find_first_of(delimiters, pos + 1);
        }
    }
}

inline string upper(string str) {
    transform(begin(str), end(str), begin(str), ::toupper);
    return str;
}

inline string lower(string str) {
    transform(begin(str), end(str), begin(str), ::tolower);
    return str;
}

inline int die(std::string str) {
    std::cout << str << std::endl;
    exit(1);
}

inline std::string env(std::string str) {
    char *path = getenv(str.c_str());
    std::string env_path = (path ? path : throw "Please set the 'NSP' environment variable to the path of the library of nsp!");
    return env_path;
}

} /// namespace jian

#endif





