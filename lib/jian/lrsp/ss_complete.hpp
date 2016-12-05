#pragma once

#include <array>
#include <list>
#include <string>
#include "../utils/traits.hpp"

BEGIN_JN

namespace lrsp {

using pair_t = std::array<int, 2>;
using pairs_t = std::list<pair_t>;
using seq_t = std::string;
using ss_t = std::string;

S ss_complete(S seq, S ss);
void ss_complete(pairs_t & pairs, const seq_t & seq);

} // namespace lrsp
END_JN

