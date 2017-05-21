#pragma once

#include "tsdna_module.hpp"

BEGIN_JN

namespace tsdna {

class TLoop : public TModule {
public:
    TLoop(const Tuple &, const Tuple &);
    virtual S type() const;
};

} // namespace tsdna

END_JN


