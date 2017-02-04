#pragma once

#include "jian/utils/log.hpp"
#include "Job.hpp"
#include "Gradient.hpp"

BEGIN_JN
namespace dg {

class CG : public Gradient {
public:
    CG() = default;
    CG(const CG &) = default;
    CG &operator =(const CG &) = default;

    void cg() {
        Mat d_o(len, 3), d_n(len, 3);

//        LOG << "Start CG..." <<  std::endl;
//        LOG << "step" << ' ' << "factor" << ' ' << "energy" << ' ' << "chirality" << std::endl;

        double a = 0.0001, beta = 0;
        for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) d_o(i, j) = 0;
        gradient();
        double oldE = _dist_en + _dih_en, oldG2 = g2;
        int upd = 0;
        while (upd < 500) {
            double tempG2 = g2, tempE = _dist_en + _dih_en;

            if (upd != 0) beta = g2 / oldG2;
            if (upd % 10 == 0) beta = 0;
            for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) {
                d_n(i, j) = -g(i, j) + beta * d_o(i, j);
                c(i, j) += a * d_n(i, j);
            }

            gradient();

            if (_dist_en + _dih_en >= tempE || g2 >= tempG2) {
                for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) c(i, j) -= a * d_n(i, j);
                gradient();
                a *= 0.5;
                if (a < 1.e-12) break;
            } else {
                upd++;
                a *= 2;
                for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) d_o(i, j) = d_n(i, j);
                oldG2 = tempG2;
                oldE = tempE;
//                LOG << upd << ' ' << a << ' ' << _dist_en << ' ' << _dih_en << std::endl;
                if (_dist_en + _dih_en < 1.e-12 || g2 < 1.e-12) break;
            }
        }
//        LOG << "Finish CG.\n"<< "Energy: "<< total_dist_energy(c)<< "(dist) "<< total_dih_energy()<< "(chir)\n" << std::endl;
    }

};

} // namespace dg
END_JN

