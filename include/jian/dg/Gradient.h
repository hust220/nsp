#ifndef JIAN_DG_GRADIENT_H
#define JIAN_DG_GRADIENT_H

#include "Job.h"
#include "../geom.h"

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
        Mat C = Mat::Zero(3 * len, 3);
        double err = 1.e-6;

        // dist energy
        _dist_en = 0;
        for (int i = 0; i < len; i++) for (int j = i + 1; j < len; j++) {
            if (i == j) continue;
            _dist_en += atom_pair_energy(c, i, j);
            for (int t = 0; t < 3; t++) {
                c(i, t) -= err; C(3 * i + t, 0) += atom_pair_energy(c, i, j); c(i, t) += err;
                c(i, t) += err; C(3 * i + t, 2) += atom_pair_energy(c, i, j); c(i, t) -= err;
                c(j, t) -= err; C(3 * j + t, 0) += atom_pair_energy(c, i, j); c(j, t) += err;
                c(j, t) += err; C(3 * j + t, 2) += atom_pair_energy(c, i, j); c(j, t) -= err;
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

    double total_energy(const Mat &coord) {
        return total_dist_energy(coord) + total_dih_energy();
    }

    double total_dist_energy(const Mat &coord) {
        double en = 0; for (int i = 0; i < coord.rows(); i++) for (int j = i + 1; j < coord.rows(); j++) {
            en += atom_pair_energy(coord, i, j);
        } return en;
    }

    double total_dih_energy() {
        double en = 0; for (auto &&entry : _dih_bound) en += dih_row_energy(entry.first); return en;
    }

    double atom_energy(const Mat &coord, int n) {
        return atom_dist_energy(coord, n) + atom_dih_energy(n);
    }

    double atom_dist_energy(const Mat &coord, int n) {
        double sum = 0; for (int i = 0; i < coord.rows(); i++) sum += atom_pair_energy(coord, i, n); return sum;
    }

    double atom_pair_energy(const Mat &coord, int i, int j) {
        if (i == j) return 0;
        double dist =  0; for (int k = 0; k < 3; k++) dist += square(coord(i, k)-coord(j, k)); dist = std::sqrt(dist);
        return (i < j ? relative_dist_energy(dist, (bound(j, i)), (bound(i, j))) : 
                        relative_dist_energy(dist, (bound(i, j)), (bound(j, i))));
    }

    double relative_dist_energy(double dist, double l, double u) {
        if (dist > u) return square(dist - u); else if (dist < l) return square(l - dist); else return 0;
    }

    double atom_dih_energy(int n) {
        return fold([&](double sum, const auto &ls){return sum + this->dih_row_energy(ls);}, 0, _involved_dihs[n]);
    }

    template<typename LS>
    double dih_row_energy(const LS &ls) {
        double dih = geom::dihedral(c.row(ls[0]), c.row(ls[1]), c.row(ls[2]), c.row(ls[3]));
        auto dih_bound = _dih_bound[ls];
        return relative_dih_energy(dih, (dih_bound.first), (dih_bound.second));
    }

    double relative_dih_energy(double en, double l, double u) {
        std::function<double(double)> foo = [&foo](auto n){
            if (n >= 3.14) return foo(n - 3.14); 
            else if (n < -3.14) return foo(n + 3.14); 
            else return n;};
        if (l < u) {
            if (l < en and en < u) return 0;
            else return std::min(foo(fabs(l - en)), foo(fabs(u - en))) / (2 * PI);
        } else if (l > u) {
            if (l < en and en < u) return std::min(foo(en - l), foo(u - en)) / (2 * PI);
            else return 0;
        } else {
            return foo(fabs(en - l)) / (2 * PI);
        }
    }

};

} // namespace dg
} // namespace jian

#endif

