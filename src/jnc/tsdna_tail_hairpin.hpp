#pragma once

#include "tsdna_module.hpp"

namespace jian {

namespace tsdna {

class TTailHairpin : public TModule {
public:
    TTailHairpin(const Tuple &, const Tuple &);
    virtual S type() const;
};

} // namespace tsdna

}


