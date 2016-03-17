#ifndef JIAN_ETL_EXISTS
#define JIAN_ETL_EXISTS

#include "std.hpp"

#define JN_EXISTS_DATA1(l)\
    BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(l),\
        std::vector<TYPEOF(PP_TUPLE_ELEM(0, l))>({PP_VARIADIC(l)});,\
        l)

#define JN_EXISTS_DATA(state)\
    JN_EXISTS_DATA1(PP_TUPLE_ELEM(PP_TUPLE_ELEM(0, state), PP_TUPLE_ELEM(2, state)))

#define JN_EXISTS_PRED1(r, state)\
    BOOST_PP_NOT_EQUAL(PP_TUPLE_ELEM(0, state), PP_TUPLE_ELEM(1, state))

#define JN_EXISTS_OP1(r, state)\
    (BOOST_PP_INC(PP_TUPLE_ELEM(0, state)), PP_TUPLE_ELEM(1, state), PP_TUPLE_ELEM(2, state))

#define JN_EXISTS_MACRO1(r, state)\
    auto && BOOST_PP_CAT(l, PP_TUPLE_ELEM(0, state)) = JN_EXISTS_DATA(state);

#define JN_EXISTS_PRED(r, state)\
    BOOST_PP_NOT_EQUAL(PP_TUPLE_ELEM(0, state), PP_TUPLE_ELEM(1, state))

#define JN_EXISTS_OP(r, state)\
    (BOOST_PP_INC(PP_TUPLE_ELEM(0, state)), PP_TUPLE_ELEM(1, state))

#define JN_EXISTS_MACRO(r, state)\
        auto && BOOST_PP_CAT(_, BOOST_PP_INC(PP_TUPLE_ELEM(0, state))) = \
            *std::next(std::begin(BOOST_PP_CAT(l, PP_TUPLE_ELEM(0, state))), i);

#define JN_EXISTS_SIZE(l)\
    BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(l),\
        PP_TUPLE_SIZE(l),\
        l.size())

#define JN_EXISTS(f, ...) ({\
    bool r = false;\
    BOOST_PP_FOR((0, PP_VARIADIC_SIZE(__VA_ARGS__), (__VA_ARGS__)), JN_EXISTS_PRED1, JN_EXISTS_OP1, JN_EXISTS_MACRO1)\
    for (int i = 0; i < JN_EXISTS_SIZE(PP_VARIADIC_ELEM(0, __VA_ARGS__)); i++) {\
        BOOST_PP_FOR((0, PP_VARIADIC_SIZE(__VA_ARGS__)), JN_EXISTS_PRED, JN_EXISTS_OP, JN_EXISTS_MACRO)\
        if (f) {r = true; break;}\
    }\
    r;\
})

#define EXISTS(f, l) ({\
    bool r = false;\
    BOOST_PP_IF(BOOST_PP_IS_BEGIN_PARENS(l),\
    for (auto &&_1 : {PP_VARIADIC(l)}) {,\
    for (auto &&_1 : l) {)\
        if (f) {r = true; break;}\
    }\
    r;\
})

#endif

