#pragma once

#include "TModule.hpp"

namespace jian {
namespace nuc3d {
namespace triple {

class TTailHairpin : public TModule {
public:
    TTailHairpin(const Tuple &, const Tuple &);
    virtual std::string type() const;
};

} // namespace triple
}
} // namespace jian


