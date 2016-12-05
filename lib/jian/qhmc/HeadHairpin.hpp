#pragma once

#include "Module.hpp"

BEGIN_JN
namespace qhmc {

class HeadHairpin : public Module {
public:
    HeadHairpin(const Tuple &, const Tuple &, int);
    virtual S type() const;
};

} // namespace quadruple
}


