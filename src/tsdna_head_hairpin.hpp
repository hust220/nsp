#pragma once

#include "tsdna_module.hpp"

namespace jian {

namespace tsdna {

    class THeadHairpin : public TModule {
        public:
            THeadHairpin(const Tuple &, const Tuple &);
            virtual S type() const;
    };

} // namespace triple

}


