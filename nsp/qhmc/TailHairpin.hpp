#pragma once

#include "Module.hpp"

BEGIN_JN
namespace qhmc {

class TailHairpin : public Module {
public:
    TailHairpin(const Tuple &, const Tuple &, int);
    virtual S type() const;
};

} // namespace quadruple
}


