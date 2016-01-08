#ifndef JIAN_ETL_LS
#define JIAN_ETL_LS

#include "core.h"

namespace jian {
namespace etl {

template<typename LS>
inline int len(LS &&ls) {
    return std::distance(std::begin(ls), std::end(ls));
}

template<typename LS, std::enable_if_t<is_same_template<LS, std::list>::value, int> = 42>
inline value_type_t<LS> &ref(LS &&ls, int n) {
    if (n < 0) return *std::next(ls.end(), n);
    else return *std::next(ls.begin(), n);
}

template<typename LS, std::enable_if_t<!is_same_template<LS, std::list>::value, int> = 42>
inline value_type_t<LS> &ref(LS &&ls, int n) {
    if (n < 0) return ls[len(ls) + n];
    else return ls[n];
}

} // namespace etl
} // namespace jian

#endif





