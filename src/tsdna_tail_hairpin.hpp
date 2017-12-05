#pragma once

#include "tsdna_module.hpp"

BEGIN_JN

namespace tsdna {

class TTailHairpin : public TModule {
public:
    TTailHairpin(const Tuple &, const Tuple &);
    virtual S type() const;
};

} // namespace tsdna

END_JN


