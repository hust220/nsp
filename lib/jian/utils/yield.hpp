#pragma once

#include "traits.hpp"

#define itBegin switch(__state__) { case 0:

#define itYield(x) do { __state__=__LINE__; return x; \
    case __LINE__:; } while (0)

#define itReturn(x) do { \
    case __LINE__: \
    __state__ = __LINE__; \
    return x; \
} while (0)

#define itEnd }

#define for_(a, b) do { auto __jn_it_##a = b; \
    for (auto a = __jn_it_##a(); a != NULL; a = __jn_it_##a())

#define end_ } while(0)

BEGIN_JN

struct itBase {
    int __state__;
    itBase(int n = 0) : __state__(n) {}
};

template<typename _Ls>
Int it_size(const _Ls &ls) {
    Int n = 0; for_(i, ls) n++; end_;
    return n;
}

// map of stl
template<typename _Out, typename _Fn, typename _In>
_Out map_j(_Fn && f, _In && in) {
    _Out out;
    for_ (i, in) {
        out.push_back(*i);
    } end_;
    return out;
}

END_JN

