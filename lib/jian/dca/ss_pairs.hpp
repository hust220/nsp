#pragma once

#include <string>
#include <list>
#include <array>

namespace jian {
namespace dca {

using pair_t = std::array<int, 2>;
using pairs_t = std::list<pair_t>;
using seq_t = std::string;
using ss_t = std::string;

pairs_t pairs_from_file(const std::string &file_name, int size = -1);
pairs_t pairs_from_ss(const ss_t &ss);
ss_t pairs_to_ss(const pairs_t &pairs, int size);

void pairs_sort(pairs_t &pairs);
void print_pairs(const pairs_t & pairs);

} // namespace lrsp
} // namespace jian

