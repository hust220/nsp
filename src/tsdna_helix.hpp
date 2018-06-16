#pragma once

#include "tsdna_module.hpp"

namespace jian {

namespace tsdna {

    class THelix : public TModule {
        public:
            THelix(const Tuple &, const Tuple &);
            virtual S type() const;
    };

} // namespace tsdna

}


