#pragma once

#include <array>
#include <list>
#include <string>

namespace jian {
namespace lrsp {

using pair_t = std::array<int, 2>;
using pairs_t = std::list<pair_t>;
using seq_t = std::string;
using ss_t = std::string;

std::string ss_complete(std::string seq, std::string ss);
void ss_complete(pairs_t & pairs, const seq_t & seq);

} // namespace lrsp
} // namespace jian

