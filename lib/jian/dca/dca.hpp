#pragma once

#include <list>
#include <string>
#include "Mat4.hpp"

namespace jian {
namespace dca {

struct mi_t {int n1, n2; double mi;};
struct di_t {int n1, n2; double di;};
using mis_t = std::list<mi_t>;
using dis_t = std::list<di_t>;
struct result_t {mis_t mis; dis_t dis;};

void analyze(std::string file, int n, result_t &rt);

void print_mis(const mis_t &mis);

void print_dis(const dis_t &dis);

} // namespace dca
} // namespace jian


