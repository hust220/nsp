#ifndef JIAN_FPL_CORE
#define JIAN_FPL_CORE

#include "../util/std.h"

namespace jian {
namespace fpl {

template<typename Fn>
auto let(Fn &&f) {
    return f();
}

} // namespace fpl
} // namespace jian

#endif





