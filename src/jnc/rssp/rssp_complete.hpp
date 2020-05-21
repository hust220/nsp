#pragma once

#include <array>
#include <list>
#include <string>
#include "jian.hpp"

namespace jian {

namespace lrsp {

using pair_t = Array<int, 2>;
using pairs_t = Deque<pair_t>;
using seq_t = Str;
using ss_t = Str;

Str ss_complete(Str seq, Str ss);
void ss_complete(pairs_t & pairs, const seq_t & seq);

} // namespace lrsp
}

