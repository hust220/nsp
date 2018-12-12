#pragma once

#include "tsdna_module.hpp"

namespace jian {

namespace tsdna {

class TLoop : public TModule {
public:
    TLoop(const Tuple &, const Tuple &);
    virtual S type() const;
};

} // namespace tsdna

}


