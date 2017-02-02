#pragma once

#include "Module.hpp"

BEGIN_JN
namespace qhmc {

class Loop : public Module {
public:
    Loop(const Tuple &, const Tuple &, int);
    virtual S type() const;
};

} // namespace quadruple
}


