#pragma once

//#include <boost/preprocessor.hpp>

// Macro
#define TYPEOF decltype
//#define CLASS_NAME(T) #T
#define PP_CALL(f, m) f m
#define PP_STRING(T) #T
#define PP_STRING1(T) PP_STRING(T)
#define PP_STRING2(T) PP_STRING1(T)
#define PP_STRING3(T) PP_STRING2(T)
#define PP_STRING4(T) PP_STRING3(T)
#define PP_STRING5(T) PP_STRING4(T)
#define PP_CAT(a, b) a##b
#define PP_CAT1(a, b) PP_CAT(a, b)
#define PP_CAT2(a, b) PP_CAT1(a, b)
#define PP_CAT3(a, b) PP_CAT2(a, b)
#define PP_CAT4(a, b) PP_CAT3(a, b)
#define PP_CAT5(a, b) PP_CAT4(a, b)
#define PP_IDENTITY(...) __VA_ARGS__
#define PP_VARIADIC(l) PP_IDENTITY l
#define PP_EMPTY(...) 
#define PP_VARIADIC_ELEM_0(e1, ...) e1
#define PP_VARIADIC_ELEM_1(e1, e2, ...) e2
#define PP_VARIADIC_ELEM_2(e1, e2, e3, ...) e3
#define PP_VARIADIC_ELEM_3(e1, e2, e3, e4, ...) e4
#define PP_VARIADIC_ELEM_4(e1, e2, e3, e4, e5, ...) e5
#define PP_VARIADIC_ELEM_5(e1, e2, e3, e4, e5, e6, ...) e6
#define PP_VARIADIC_ELEM(i, ...) PP_CAT2(PP_VARIADIC_ELEM_, i)(__VA_ARGS__, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#define PP_TUPLE_ELEM(i, n) PP_VARIADIC_ELEM(i, PP_VARIADIC(n))
#define PP_VARIADIC_SIZE_HELPER(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15, e16, e17, e18, e19, e20, e21, e22, e23, e24, e25, e26, e27, e28, e29, e30, size, ...) size
#define PP_VARIADIC_SIZE(...) PP_CALL(PP_VARIADIC_SIZE_HELPER, (__VA_ARGS__, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1))
#define PP_TUPLE_SIZE(i) PP_VARIADIC_SIZE i

//// IF
//#define JN_IF(...) do{\
//    BOOST_PP_IF(BOOST_PP_EQUAL(PP_VARIADIC_SIZE(__VA_ARGS__), 2),\
//        if (PP_VARIADIC_ELEM(0, __VA_ARGS__)) {PP_VARIADIC_ELEM(1, __VA_ARGS__);},\
//        if (PP_VARIADIC_ELEM(0, __VA_ARGS__)) {PP_VARIADIC_ELEM(1, __VA_ARGS__);}\
//        else {PP_VARIADIC_ELEM(2, __VA_ARGS__);})\
//}while(0)
//#define IF JN_IF
//
//// SWITCH
//#define JN_SWITCH_HELPER(r, data, elem)\
//    else if (data == PP_TUPLE_ELEM(0, elem)) {PP_TUPLE_ELEM(1, elem);}
//
//#define JN_SWITCH(i, ...)\
//    if (i == PP_TUPLE_ELEM(0, BOOST_PP_LIST_FIRST(BOOST_PP_VARIADIC_TO_LIST(__VA_ARGS__))))\
//        {PP_TUPLE_ELEM(1, BOOST_PP_LIST_FIRST(BOOST_PP_VARIADIC_TO_LIST(__VA_ARGS__)));}\
//    BOOST_PP_LIST_FOR_EACH(JN_SWITCH_HELPER, i, BOOST_PP_LIST_REST(BOOST_PP_VARIADIC_TO_LIST(__VA_ARGS__)))
//
//#define SWITCH JN_SWITCH
//
//// COND
//#define JN_COND_HELPER(r, data, elem)\
//    else if (PP_TUPLE_ELEM(0, elem)) {PP_TUPLE_ELEM(1, elem);}
//
//#define JN_COND(...)\
//    if (PP_TUPLE_ELEM(0, PP_VARIADIC_ELEM(0, __VA_ARGS__)))\
//        {PP_TUPLE_ELEM(1, PP_VARIADIC_ELEM(1, __VA_ARGS__));}\
//    BOOST_PP_LIST_FOR_EACH(JN_COND_HELPER, , BOOST_PP_LIST_REST(BOOST_PP_VARIADIC_TO_LIST(__VA_ARGS__)))
//
//#define COND JN_COND

// HAS
//#define JN_HAS_HELPER(r, data, elem)\
//    || data == elem
//
//#define JN_HAS(l, i) (\
//    BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(l),\
//        i == PP_TUPLE_ELEM(0, l) BOOST_PP_LIST_FOR_EACH(JN_HAS_HELPER, i, BOOST_PP_LIST_REST(BOOST_PP_TUPLE_TO_LIST(l))),\
//        std::count(std::begin(l), std::end(l), i))\
//)
//
//#define HAS JN_HAS
//
// FOR
//#define JN_FOR_HELPER(i, beg, end, step, c) do{\
//    int i = beg; \
//    for (; i < end; i += step) {c;} \
//    i;\
//}while(0)
//
//#define JN_FOR(i, c) \
//    BOOST_PP_IF(BOOST_PP_EQUAL(PP_TUPLE_SIZE(i), 4),\
//        JN_FOR_HELPER(PP_TUPLE_ELEM(0, i), PP_TUPLE_ELEM(1, i),\
//                        PP_TUPLE_ELEM(2, i), PP_TUPLE_ELEM(3, i), c),\
//        BOOST_PP_IF(BOOST_PP_EQUAL(PP_TUPLE_SIZE(i), 2),\
//            JN_FOR_HELPER(PP_TUPLE_ELEM(0, i), 0, PP_TUPLE_ELEM(1, i), 1, c),\
//            BOOST_PP_IF(BOOST_PP_EQUAL(PP_TUPLE_SIZE(i), 3),\
//                JN_FOR_HELPER(PP_TUPLE_ELEM(0, i), PP_TUPLE_ELEM(1, i), PP_TUPLE_ELEM(2, i), 1, c),\
//                BOOST_PP_EMPTY())))
//
//#define FOR JN_FOR
//
// EACH
//#define JN_EACH_HELPER(i, l, c) do{\
//    BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(l),\
//        for (auto &&i : {PP_VARIADIC(l)}) {c;},\
//        for (auto &&i : l) {c;})\
//}while(0)
//
//#define JN_INDEX_EACH_HELPER(i, n, l, c) do{\
//    int n = 0; \
//    BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(l),\
//        for (auto &&i : {PP_VARIADIC(l)}) {c; n++;},\
//        for (auto &&i : l) {c; n++;}) \
//    n;\
//}while(0)
//
//#define JN_EACH(i, l, c) \
//    BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(i),\
//        JN_INDEX_EACH_HELPER(PP_TUPLE_ELEM(0, i), PP_TUPLE_ELEM(1, i), l, c),\
//        JN_EACH_HELPER(i, l, c))
//
//#define EACH JN_EACH
//
//// LOOP
//#define JN_LOOP(i, n, c) ({for (int i = 0; i < n; i++) {c;}})
//#define LOOP JN_LOOP

// FOLD
//#define JN_FOLD(f, n, l) do{\
//        BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(n), \
//            PP_TUPLE_ELEM(0, n) _1 = PP_TUPLE_ELEM(1, n);, \
//            auto _1 = n;)\
//        BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(l), \
//            for (auto &&_2 : {PP_VARIADIC(l)}) {, \
//            for (auto &&_2 : l) {) \
//            _1 = ({f;}); \
//        }\
//        _1;\
//}while(0)
//
//#define FOLD JN_FOLD
//
//#define JN_EXISTS_DATA1(l)\
//    BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(l),\
//        std::vector<TYPEOF(PP_TUPLE_ELEM(0, l))>({PP_VARIADIC(l)});,\
//        l)
//
//#define JN_EXISTS_DATA(state)\
//    JN_EXISTS_DATA1(PP_TUPLE_ELEM(PP_TUPLE_ELEM(0, state), PP_TUPLE_ELEM(2, state)))
//
//#define JN_EXISTS_PRED1(r, state)\
//    BOOST_PP_NOT_EQUAL(PP_TUPLE_ELEM(0, state), PP_TUPLE_ELEM(1, state))
//
//#define JN_EXISTS_OP1(r, state)\
//    (BOOST_PP_INC(PP_TUPLE_ELEM(0, state)), PP_TUPLE_ELEM(1, state), PP_TUPLE_ELEM(2, state))
//
//#define JN_EXISTS_MACRO1(r, state)\
//    auto && BOOST_PP_CAT(l, PP_TUPLE_ELEM(0, state)) = JN_EXISTS_DATA(state);
//
//#define JN_EXISTS_PRED(r, state)\
//    BOOST_PP_NOT_EQUAL(PP_TUPLE_ELEM(0, state), PP_TUPLE_ELEM(1, state))
//
//#define JN_EXISTS_OP(r, state)\
//    (BOOST_PP_INC(PP_TUPLE_ELEM(0, state)), PP_TUPLE_ELEM(1, state))
//
//#define JN_EXISTS_MACRO(r, state)\
//        auto && BOOST_PP_CAT(_, BOOST_PP_INC(PP_TUPLE_ELEM(0, state))) = \
//            *std::next(std::begin(BOOST_PP_CAT(l, PP_TUPLE_ELEM(0, state))), i);
//
//#define JN_EXISTS_SIZE(l)\
//    BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(l),\
//        PP_TUPLE_SIZE(l),\
//        l.size())
//
//#define JN_EXISTS(f, ...) do{\
//    bool r = false;\
//    BOOST_PP_FOR((0, PP_VARIADIC_SIZE(__VA_ARGS__), (__VA_ARGS__)), JN_EXISTS_PRED1, JN_EXISTS_OP1, JN_EXISTS_MACRO1)\
//    for (int i = 0; i < JN_EXISTS_SIZE(PP_VARIADIC_ELEM(0, __VA_ARGS__)); i++) {\
//        BOOST_PP_FOR((0, PP_VARIADIC_SIZE(__VA_ARGS__)), JN_EXISTS_PRED, JN_EXISTS_OP, JN_EXISTS_MACRO)\
//        if (f) {r = true; break;}\
//    }\
//    r;\
//}while(0)
//
//#define EXISTS(f, l) do{\
//    bool r = false;\
//    BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(l),\
//    for (auto &&_1 : {PP_VARIADIC(l)}) {,\
//    for (auto &&_1 : l) {)\
//        if (f) {r = true; break;}\
//    }\
//    r;\
//}while(0)
