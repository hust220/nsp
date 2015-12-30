#ifndef JIAN_DG_GRADIENT_H
#define JIAN_DG_GRADIENT_H

#include "Job.h"
#include "../geom/geometry.h"

namespace jian {
namespace dg {

class Gradient : public virtual Job {
public:
    std::vector<std::list<DihEntry>> _involved_dihs;

    Gradient() {
        init_dihs();
    }

    Gradient(const Gradient &) = default;
    Gradient &operator =(const Gradient &) = default;

    void gradient() {
        MatrixXf C = MatrixXf::Zero(3 * len, 3);
        double err = 1.e-6;

        // dist energy
        _dist_en = 0;
        for (int i = 0; i < len; i++) for (int j = i + 1; j < len; j++) {
            if (i == j) continue;
            _dist_en += atom_pair_energy(c, i, j);
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
        _dih_en = 0;
        for (auto && pair : _dih_bound) {
            _dih_en += dih_row_energy(pair.first);
            for (auto && index : pair.first) for (int t = 0; t < 3; t++) for (auto l : {0, 2}) {
                c(index, t) += (l - 1) * err; C(index * 3 + t, l) += dih_row_energy(pair.first);
                c(index, t) -= (l - 1) * err;
            }
        }

        g.resize(len, 3); g2 = 0;
        for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) {
            g(i, j) = (C(3 * i + j, 2) - C(3 * i + j, 0)) / (2 * err);
            g2 += g(i, j) * g(i, j);
        }
    }

    void init_dihs() {
        _involved_dihs.resize(c.rows());
        for (auto && pair : _dih_bound) for (auto i : pair.first) _involved_dihs[i].push_back(pair.first);
    }

    double total_energy(const MatrixXf &coord) {
        return total_dist_energy(coord) + total_dih_energy();
    }

    double total_dist_energy(const MatrixXf &coord) {
        double en = 0, u2, l2, dist2;
        for (int i = 0; i < coord.rows(); i++) {
            for (int j = i + 1; j < coord.rows(); j++) {
                dist2 = fold([&](double n, int k)->double{return n + square(coord(i, k) - coord(j, k));}, 0, {0, 3});
                (u2 = square(bound(i, j)), l2 = square(bound(j, i)));
                if (dist2 > u2) en += (dist2 - u2) / u2;
                else if (dist2 < l2) en += (l2 - dist2) / dist2;
            }
        }
        return en;
    }

    double total_dih_energy() {
        return fold([&](double sum, const auto &entry){return sum + this->dih_row_energy(entry.first);}, 0, _dih_bound);
    }

    double atom_energy(const MatrixXf &coord, int n) {
        return atom_dist_energy(coord, n) + atom_dih_energy(n);
    }

    double atom_dist_energy(const MatrixXf &coord, int n) {
        return fold([&](double sum, auto i){return sum + this->atom_pair_energy(coord, i, n);}, 0, {0, coord.rows()});
    }

    double atom_pair_energy(const MatrixXf &coord, int i, int j) {
        if (i == j) return 0;
        double dist2 = fold([&](double sum, int k){return sum + square(coord(i, k) - coord(j, k));}, 0, {0, 3});
        return (i < j ? relative_energy(dist2, square(bound(j, i)), square(bound(i, j))) : 
                        relative_energy(dist2, square(bound(i, j)), square(bound(j, i))));
    }

    double atom_dih_energy(int n) {
        return fold([&](double sum, const auto &ls){return sum + this->dih_row_energy(ls);}, 0, _involved_dihs[n]);
    }

    template<typename LS>
    double dih_row_energy(const LS &ls) {
        double dih2 = geometry::chirality(c.row(ls[0]), c.row(ls[1]), c.row(ls[2]), c.row(ls[3]));
        auto dih_bound = _dih_bound[ls];
        return relative_energy(dih2, square(dih_bound.first), square(dih_bound.second));
    }

    double relative_energy(double en2, double l2, double u2) {
        if (en2 > u2) return (en2 - u2) / u2; else if (en2 < l2) return (l2 - en2) / en2; else return 0;
    }

};

} // namespace dg
} // namespace jian

#endif


