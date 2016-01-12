#ifndef JIAN_DG_CG
#define JIAN_DG_CG

#include "Sys.h"

namespace jian {
namespace dg {

class CG {
public:
    template<typename Sys> operator ()(Sys &&sys) {
        cg(sys);
    }

    template<typename Sys> void cg(Sys &&sys) {
        MatrixXf old_direction = MatrixXf::Zero(len, 3), new_direction(len, 3);
        auto old_state = sys.state(State::EN | State::G | State::G2);
        double factor = 0.001, beta = 0;
        int upd = 0;

        log("Start CG...\n", 
            "Energy: ", old_state.energy, "\n");

        while (upd < 500) {
            if (upd != 0) beta = new_state.g2 / old_state.g2;
            if (upd % 10 == 0) beta = 0;
            for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) {
                new_direction(i, j) = -old_state.g(i, j) + beta * old_direction(i, j);
                sys.c(i, j) += factor * new_direction(i, j);
            }

            auto new_state = sys.state(State::EN | State::G | State::G2);

            if (new_state.energy >= old_state.energy or new_state.g2 >= old_state.g2) {
                for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) {
                    sys.c(i, j) -= factor * new_direction(i, j);
                }
                factor *= 0.5;
                if (factor < 1.e-12) break;
            } else {
                upd++;
                factor *= 2;
                old_direction = new_direction;
                old_state = std::move(new_state);
                log(upd, ' ', factor, ' ', old_state.energy, '\n');
                if (old_state.energy < 1.e-12 and old_state.g2 < 1.e-12) break;
            }
        }
        log("Finish CG.\n",
            "Energy: ", old_state.energy, "\n\n");
    }

};

} // namespace dg
} // namespace jian

#endif


