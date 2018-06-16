#pragma once

#include "g4dna_module.hpp"

namespace jian {
namespace qhmc {

class Helix : public Module {
public:
    Helix(const Tuple &, const Tuple &, int);
    virtual S type() const;
};

} // namespace quadruple
}


