#pragma once

#include <string>
#include <list>
#include <array>
#include  "../utils/string.hpp"

namespace jian {
namespace dca {

struct tuple_t {
    int a, b;
    double c;
};
using tuples_t = std::list<tuple_t>;
using pair_t = std::array<int, 2>;
using pairs_t = std::list<pair_t>;
using seq_t = str_t;
using ss_t = str_t;

tuples_t tuples_from_file(const str_t &file_name, int size);
pairs_t pairs_from_file(const str_t &file_name, int size);
pairs_t pairs_from_ss(const ss_t &ss);
ss_t pairs_to_ss(const pairs_t &pairs, int size);

void pairs_sort(pairs_t &pairs);
void print_pairs(const pairs_t & pairs);
void print_tuples(const tuples_t & tuples);

} // namespace lrsp
} // namespace jian

