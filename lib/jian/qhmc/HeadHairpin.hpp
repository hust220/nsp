#pragma once

#include "Module.hpp"

namespace jian {
namespace qhmc {

class HeadHairpin : public Module {
public:
    HeadHairpin(const Tuple &, const Tuple &, int);
    virtual std::string type() const;
};

} // namespace quadruple
}


