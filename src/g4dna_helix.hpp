#pragma once

#include "g4dna_module.hpp"

BEGIN_JN
namespace qhmc {

class Helix : public Module {
public:
    Helix(const Tuple &, const Tuple &, int);
    virtual S type() const;
};

} // namespace quadruple
}


