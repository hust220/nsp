#ifndef JIAN_MLIB_H
#define JIAN_MLIB_H

#include "std.h"

using namespace std;

namespace jian {

template <typename T, size_t N>
char (&_ArraySizeHelper(T (&array)[N]))[N];

template <typename Array>
std::size_t count(Array array) {
    return sizeof _ArraySizeHelper(array);
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

void tokenize(const string &, vector<string> &, const string & = " ");
void tokenize(const string &, vector<string> &, const string &, const string &temp);
string upper(string);
string lower(string);
int die(string str);
std::string env(std::string);

} /// namespace jian

#endif





