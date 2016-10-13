#pragma once

#include "Module.hpp"

namespace jian {
namespace qhmc {

class Helix : public Module {
public:
    Helix(const Tuple &, const Tuple &);
    virtual std::string type() const;
};

} // namespace quadruple
}


