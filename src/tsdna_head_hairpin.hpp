#pragma once

#include "tsdna_module.hpp"

BEGIN_JN

namespace tsdna {

    class THeadHairpin : public TModule {
        public:
            THeadHairpin(const Tuple &, const Tuple &);
            virtual S type() const;
    };

} // namespace triple

END_JN


