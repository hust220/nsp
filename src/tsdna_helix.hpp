#pragma once

#include "tsdna_module.hpp"

BEGIN_JN

namespace tsdna {

    class THelix : public TModule {
        public:
            THelix(const Tuple &, const Tuple &);
            virtual S type() const;
    };

} // namespace tsdna

END_JN


