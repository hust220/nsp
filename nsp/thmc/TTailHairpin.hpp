#pragma once

#include "TModule.hpp"

BEGIN_JN
namespace nuc3d {
namespace triple {

class TTailHairpin : public TModule {
public:
    TTailHairpin(const Tuple &, const Tuple &);
    virtual S type() const;
};

} // namespace triple
}
END_JN


