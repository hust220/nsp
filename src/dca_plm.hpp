#pragma once

#include "dca.hpp"

namespace jian {

namespace dca {

    class PlmDca : public Dca {
        public:
            Matf C, eij;

            PlmDca();

            PlmDca(S mol_type, float pw);

            virtual void calculate_eij();

            virtual float cal_di(int i, int j);
    };

} // namespace dca

}


