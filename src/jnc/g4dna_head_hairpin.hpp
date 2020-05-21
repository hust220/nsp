#pragma once

#include "g4dna_module.hpp"

namespace jian {
namespace qhmc {

class HeadHairpin : public Module {
public:
    HeadHairpin(const Tuple &, const Tuple &, int);
    virtual S type() const;
};

} // namespace quadruple
}


