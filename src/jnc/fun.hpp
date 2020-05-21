#pragma once

//#define SIZE(...) BOOST_PP_TUPLE_SIZE((__VA_ARGS__))
//
//#define _CAT(a, b) a##b
//#define CAT(a, b) _CAT(a, b)
//
//#define _STRING(a) #a
//#define STRING(a) _STRING(a)
//
//#define OP_HEAD(d, state, x) \
//    (BOOST_PP_SUB(BOOST_PP_TUPLE_ELEM(2, 0, state), 1), \
//     BOOST_PP_LIST_CONS(typename CAT(_, BOOST_PP_TUPLE_ELEM(2, 0, state)), BOOST_PP_TUPLE_ELEM(2, 1, state)))
//
//#define OP_BODY(d, state, x) \
//    (BOOST_PP_SUB(BOOST_PP_TUPLE_ELEM(2, 0, state), 1), \
//     BOOST_PP_LIST_CONS(CAT(_, BOOST_PP_TUPLE_ELEM(2, 0, state)) x, BOOST_PP_TUPLE_ELEM(2, 1, state)))
//
//#define FUN_HEAD(...) \
//    template<BOOST_PP_REMOVE_PARENS( \
//        BOOST_PP_LIST_TO_TUPLE( \
//            BOOST_PP_TUPLE_ELEM(2, 1, \
//                BOOST_PP_LIST_FOLD_RIGHT(OP_HEAD, (SIZE(__VA_ARGS__), BOOST_PP_LIST_NIL), \
//                    BOOST_PP_TUPLE_TO_LIST((__VA_ARGS__))))))>
//
//#define FUN_BODY(...) \
//    BOOST_PP_LIST_TO_TUPLE( \
//        BOOST_PP_TUPLE_ELEM(2, 1, \
//            BOOST_PP_LIST_FOLD_RIGHT(OP_BODY, (SIZE(__VA_ARGS__), BOOST_PP_LIST_NIL), \
//                BOOST_PP_TUPLE_TO_LIST((__VA_ARGS__)))))
//
//#define FUN(f, ...) \
//    FUN_HEAD(__VA_ARGS__) auto f FUN_BODY(__VA_ARGS__)
//
