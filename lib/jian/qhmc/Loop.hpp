#pragma once

#include "Module.hpp"

namespace jian {
namespace qhmc {

class Loop : public Module {
public:
    Loop(const Tuple &, const Tuple &);
    virtual std::string type() const;
};

} // namespace quadruple
}


