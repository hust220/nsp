#pragma once

#include "Dca.hpp"

BEGIN_JN

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

END_JN


