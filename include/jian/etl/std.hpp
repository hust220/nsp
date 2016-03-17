#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstring>
#include <cassert>
#include <regex>

// STL
#include <bitset>
#include <vector>
#include <array>
#include <queue>
#include <numeric>
#include <ctime>
#include <list>
#include <set>
#include <deque>
#include <unordered_map>
#include <map>
#include <memory>
#include <iterator>
#include <utility>
#include <algorithm>
#include <functional>
#include <type_traits>
using namespace std;

// boost
#include <boost/timer.hpp>
#include <boost/assert.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/preprocessor.hpp>
#include <boost/lexical_cast.hpp>

// Macro
#define MASS_H 1.00794
#define MASS_C 12.0107
#define MASS_O 15.9994
#define MASS_N 14.0067
#define MASS_P 30.973762

#define PI 3.1415927

#define TYPEOF decltype
#define CLASS_NAME(T) #T
#define PP_STRING(T) #T

#define PP_CAT(a, b) a##b
#define PP_IDENTITY(...) __VA_ARGS__
#define PP_VARIADIC(l) PP_IDENTITY l
#define PP_EMPTY(...) 
#define PP_VARIADIC_ELEM_0(e1, ...) e1
#define PP_VARIADIC_ELEM_1(e1, e2, ...) e2
#define PP_VARIADIC_ELEM_2(e1, e2, e3, ...) e3
#define PP_VARIADIC_ELEM_3(e1, e2, e3, e4...) e4
#define PP_VARIADIC_ELEM_4(e1, e2, e3, e4, e5...) e5
#define PP_VARIADIC_ELEM_5(e1, e2, e3, e4, e5, e6...) e6
#define PP_VARIADIC_ELEM(i, ...) BOOST_PP_CAT(PP_VARIADIC_ELEM_, i)(__VA_ARGS__, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#define PP_TUPLE_ELEM(i, n) PP_VARIADIC_ELEM(i, PP_VARIADIC(n))
#define PP_VARIADIC_SIZE_HELPER(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, size, ...) size
#define PP_VARIADIC_SIZE(...) PP_VARIADIC_SIZE_HELPER(__VA_ARGS__, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
#define PP_TUPLE_SIZE(i) PP_VARIADIC_SIZE i

// IF
#define JN_IF(...) ({\
    BOOST_PP_IF(BOOST_PP_EQUAL(PP_VARIADIC_SIZE(__VA_ARGS__), 2),\
        if (PP_VARIADIC_ELEM(0, __VA_ARGS__)) {PP_VARIADIC_ELEM(1, __VA_ARGS__);},\
        if (PP_VARIADIC_ELEM(0, __VA_ARGS__)) {PP_VARIADIC_ELEM(1, __VA_ARGS__);}\
        else {PP_VARIADIC_ELEM(2, __VA_ARGS__);})\
})
#define IF JN_IF

// SWITCH
#define JN_SWITCH_HELPER(r, data, elem)\
    else if (data == PP_TUPLE_ELEM(0, elem)) {PP_TUPLE_ELEM(1, elem);}

#define JN_SWITCH(i, ...)\
    if (i == PP_TUPLE_ELEM(0, BOOST_PP_LIST_FIRST(BOOST_PP_VARIADIC_TO_LIST(__VA_ARGS__))))\
        {PP_TUPLE_ELEM(1, BOOST_PP_LIST_FIRST(BOOST_PP_VARIADIC_TO_LIST(__VA_ARGS__)));}\
    BOOST_PP_LIST_FOR_EACH(JN_SWITCH_HELPER, i, BOOST_PP_LIST_REST(BOOST_PP_VARIADIC_TO_LIST(__VA_ARGS__)))

#define SWITCH JN_SWITCH

// COND
#define JN_COND_HELPER(r, data, elem)\
    else if (PP_TUPLE_ELEM(0, elem)) {PP_TUPLE_ELEM(1, elem);}

#define JN_COND(...)\
    if (PP_TUPLE_ELEM(0, PP_VARIADIC_ELEM(0, __VA_ARGS__)))\
        {PP_TUPLE_ELEM(1, PP_VARIADIC_ELEM(1, __VA_ARGS__));}\
    BOOST_PP_LIST_FOR_EACH(JN_COND_HELPER, , BOOST_PP_LIST_REST(BOOST_PP_VARIADIC_TO_LIST(__VA_ARGS__)))

#define COND JN_COND

// HAS
#define JN_HAS_HELPER(r, data, elem)\
    || data == elem

#define JN_HAS(l, i) (\
    BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(l),\
        i == PP_TUPLE_ELEM(0, l) BOOST_PP_LIST_FOR_EACH(JN_HAS_HELPER, i, BOOST_PP_LIST_REST(BOOST_PP_TUPLE_TO_LIST(l))),\
        std::count(std::begin(l), std::end(l), i))\
)

#define HAS JN_HAS

// FOR
#define JN_FOR_HELPER(i, beg, end, step, c) ({\
    int i = beg; \
    for (; i < end; i += step) {c;} \
    i;\
})

#define JN_FOR(i, c) \
    BOOST_PP_IF(BOOST_PP_EQUAL(PP_TUPLE_SIZE(i), 4),\
        JN_FOR_HELPER(PP_TUPLE_ELEM(0, i), PP_TUPLE_ELEM(1, i),\
                        PP_TUPLE_ELEM(2, i), PP_TUPLE_ELEM(3, i), c),\
        BOOST_PP_IF(BOOST_PP_EQUAL(PP_TUPLE_SIZE(i), 2),\
            JN_FOR_HELPER(PP_TUPLE_ELEM(0, i), 0, PP_TUPLE_ELEM(1, i), 1, c),\
            BOOST_PP_IF(BOOST_PP_EQUAL(PP_TUPLE_SIZE(i), 3),\
                JN_FOR_HELPER(PP_TUPLE_ELEM(0, i), PP_TUPLE_ELEM(1, i), PP_TUPLE_ELEM(2, i), 1, c),\
                BOOST_PP_EMPTY())))

#define FOR JN_FOR

// EACH
#define JN_EACH_HELPER(i, l, c) ({\
    BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(l),\
        for (auto &&i : {PP_VARIADIC(l)}) {c;},\
        for (auto &&i : l) {c;})\
})

#define JN_INDEX_EACH_HELPER(i, n, l, c) ({\
    int n = 0; \
    BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(l),\
        for (auto &&i : {PP_VARIADIC(l)}) {c; n++;},\
        for (auto &&i : l) {c; n++;}) \
    n;\
})

#define JN_EACH(i, l, c) \
    BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(i),\
        JN_INDEX_EACH_HELPER(PP_TUPLE_ELEM(0, i), PP_TUPLE_ELEM(1, i), l, c),\
        JN_EACH_HELPER(i, l, c))

#define EACH JN_EACH

// LOOP
#define JN_LOOP(i, n, c) ({for (int i = 0; i < n; i++) {c;}})
#define LOOP JN_LOOP

// FOLD
#define JN_FOLD(f, n, l) ({\
        BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(n), \
            PP_TUPLE_ELEM(0, n) _1 = PP_TUPLE_ELEM(1, n);, \
            auto _1 = n;)\
        BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(l), \
            for (auto &&_2 : {PP_VARIADIC(l)}) {, \
            for (auto &&_2 : l) {) \
            _1 = ({f;}); \
        }\
        _1;\
})

#define FOLD JN_FOLD

// JINT JSTR
#define JN_INT(n) boost::lexical_cast<int>(n)
#define JN_FLT(n) boost::lexical_cast<float>(n)
#define JN_DBL(n) boost::lexical_cast<double>(n)
#define JN_STR(n) boost::lexical_cast<std::string>(n)
#define JINT JN_INT
#define JSTR JN_STR

inline void SEE_() {}
template<typename T, typename... U> 
inline void SEE_(T &&data, U && ...datum) { std::cout << data; SEE_(datum...); }
template<typename... U> 
inline void SEELN_(U && ...datum) {SEE_(datum...); std::cout << std::endl; }
#define SEE(...) SEE_(__VA_ARGS__)
#define SEELN(...) SEELN_(__VA_ARGS__)

