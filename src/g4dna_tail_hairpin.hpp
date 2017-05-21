#pragma once

#include "g4dna_module.hpp"

BEGIN_JN
namespace qhmc {

class TailHairpin : public Module {
public:
    TailHairpin(const Tuple &, const Tuple &, int);
    virtual S type() const;
};

} // namespace quadruple
}


