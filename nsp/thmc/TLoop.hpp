#pragma once

#include "TModule.hpp"

BEGIN_JN
namespace nuc3d {
namespace triple {

class TLoop : public TModule {
public:
    TLoop(const Tuple &, const Tuple &);
    virtual S type() const;
};

} // namespace triple
}
END_JN


