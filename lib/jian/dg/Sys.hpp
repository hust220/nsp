#pragma once

BEGIN_JN
namespace dg {

class State {
public:
    const static int EN = 1;
    const static int G = 2;
    const static int G2 = 4;

    double energy;
    MatrixXf g;
    double g2;
};

class Sys {
public:
    MatrixXf c;
    Constraints constraints;

    State state(int type = 1) {
        State s;
        MatrixXf C = MatrixXf::Zero(3 * len, 3);
        double err = 1.e-6;

        double energy = 0;
        for (int i = 0; i < len; i++) for (int j = i + 1; j < len; j++) {
            if (i == j) continue;
            double center = atom_pair_energy(c, i, j);
            energy += center;
            for (int t = 0; t < 3; t++) {
                c(i, t) -= err; C(3 * i + t, 0) += atom_pair_energy(c, i, j);
                c(i, t) += 2 * err; C(3 * i + t, 2) += atom_pair_energy(c, i, j);
                c(i, t) -= err;
                c(j, t) -= err; C(3 * j + t, 0) += atom_pair_energy(c, i, j);
                c(j, t) += 2 * err; C(3 * j + t, 2) += atom_pair_energy(c, i, j);
                c(j, t) -= err;
            }
        }
        
        // chir energy
        for (int i = 0; i < chir.rows(); i++) {
            double en = chir_row_energy(i);
            energy += en;
            int k[4]; for (int j = 0; j < 4; j++) k[j] = int(chir(i, j));
            for (auto && j : k) for (int t = 0; t < 3; t++) for (auto l : {0, 2}) {
                c(j, t) += (l - 1) * err; C(j * 3 + t, l) += chir_row_energy(i);
                c(j, t) -= (l - 1) * err;
            }
        }

        MatrixXf g(len, 3);
        double g2 = 0;
        for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) {
            g(i, j) = (C(3 * i + j, 2) - C(3 * i + j, 0)) / (2 * err);
            g2 += g(i, j) * g(i, j);
        }

    }
};

} // namespace dg
END_JN

