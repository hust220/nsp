#pragma once

#include <deque>
#include <array>
#include <memory>
#include "../matrix.hpp"
#include "../utils/Factory.hpp"

namespace jian {
namespace qhmc {

using Frag = std::deque<int>;
using Frags = std::deque<Frag>;
using Tuple = std::array<int, 4>;
using Tuples = std::deque<Tuple>;
using Tree = std::deque<Tuples>;

class Module {
public:
    using cons_t = Module*(const Tuple &, const Tuple &);
    Frags d_frags;
    int d_max_len {0};
    std::unique_ptr<Mat> d_indices;

    void set_indices(int, int, int);
    virtual std::string type() const = 0;
};

#define REGISTER_QUADRUPLE_MODULE_FACTORY(name, Type) REGISTER_FACTORY(jian::qhmc::Module::cons_t, name, Type)

} // namespace qhmc
} // namespace jian


