#pragma once

#include "g4dna_module.hpp"

namespace jian {
namespace qhmc {

class Loop : public Module {
public:
    Loop(const Tuple &, const Tuple &, int);
    virtual S type() const;
};

} // namespace quadruple
}


