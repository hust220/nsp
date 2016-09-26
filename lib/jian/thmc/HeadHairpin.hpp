#pragma once

#include "Module.hpp"

namespace jian {
namespace nuc3d {
namespace triple {

class HeadHairpin : public Module {
public:
    HeadHairpin(const Tuple &, const Tuple &);
    virtual std::string type() const;
};

} // namespace triple
}
} // namespace jian


