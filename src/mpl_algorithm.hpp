#ifndef JIAN_MPL_ALGORITHM
#define JIAN_MPL_ALGORITHM

#include "list.h"

BEGIN_JN
namespace mpl {

template<template<typename...> class F, typename D, typename L>
struct fold;

template<template<typename...> class F, typename D>
struct fold<F, D, NullType> {
    typedef D type;
};

template<template<typename...> class F, typename D, typename Head, typename Tail>
struct fold<F, D, Pair<Head, Tail>> {
    typedef typename fold<F, typename F<D, Head>::type, Tail>::type type;
};

} // namespace mpl
END_JN










#endif

