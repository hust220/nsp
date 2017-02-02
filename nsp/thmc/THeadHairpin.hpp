#pragma once

#include "TModule.hpp"

BEGIN_JN
namespace nuc3d {
namespace triple {

class THeadHairpin : public TModule {
public:
    THeadHairpin(const Tuple &, const Tuple &);
    virtual S type() const;
};

} // namespace triple
}
END_JN


