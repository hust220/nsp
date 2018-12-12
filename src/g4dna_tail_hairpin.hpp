#pragma once

#include "g4dna_module.hpp"

namespace jian {
namespace qhmc {

class TailHairpin : public Module {
public:
    TailHairpin(const Tuple &, const Tuple &, int);
    virtual S type() const;
};

} // namespace quadruple
}


