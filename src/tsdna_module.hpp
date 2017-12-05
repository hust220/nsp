#pragma once

#include <deque>
#include <array>
#include <memory>
#include "matrix.hpp"
#include "factory.hpp"

BEGIN_JN

namespace tsdna {

using Frag = std::deque<int>;
using Frags = std::deque<Frag>;
using Tuple = std::array<int, 3>;
using Tuples = std::deque<Tuple>;
using Tree = std::deque<Tuples>;

class TModule {
public:
    using cons_t = TModule*(const Tuple &, const Tuple &);
    Frags d_frags;
    int d_max_len {0};
    std::unique_ptr<Mati> d_indices;

    void set_indices(int, int, int);
    virtual S type() const = 0;
};

#define REGISTER_TRIPLE_MODULE_FACTORY(name, Type) REGISTER_FACTORY(jian::tsdna::TModule::cons_t, name, Type)

} // namespace tsdna

END_JN


