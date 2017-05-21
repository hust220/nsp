#include "dg_gradient.hpp"
#include "geom.hpp"
#include <iterator>

BEGIN_JN
namespace dg {

Gradient::Gradient() { init_dihs(); }

void Gradient::gradient() {
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

//    LOG << g2 << std::endl;
}

void Gradient::init_dihs() {
    _involved_dihs.resize(c.rows());
    for (auto && pair : _dih_bound) for (auto i : pair.first) _involved_dihs[i].push_back(pair.first);
}

double Gradient::total_energy(const Mat &coord) {
    return total_dist_energy(coord) + total_dih_energy();
}

double Gradient::total_dist_energy(const Mat &coord) {
    double en = 0; for (int i = 0; i < coord.rows(); i++) for (int j = i + 1; j < coord.rows(); j++) {
        en += atom_pair_energy(coord, i, j);
    } return en;
}

double Gradient::total_dih_energy() {
    double en = 0; 
    for (auto &&entry : _dih_bound) {
        en += dih_row_energy(entry.first); 
    }
    return en;
}

double Gradient::atom_energy(const Mat &coord, int n) {
    return atom_dist_energy(coord, n) + atom_dih_energy(n);
}

double Gradient::atom_dist_energy(const Mat &coord, int n) {
    double sum = 0; for (int i = 0; i < coord.rows(); i++) sum += atom_pair_energy(coord, i, n); return sum;
}

double Gradient::atom_pair_energy(const Mat &coord, int i, int j) {
    if (i == j) return 0;
    double dist =  0; for (int k = 0; k < 3; k++) dist += square(coord(i, k)-coord(j, k)); dist = std::sqrt(dist);
    return (i < j ? relative_dist_energy(dist, (bound(j, i)), (bound(i, j))) : 
                    relative_dist_energy(dist, (bound(i, j)), (bound(j, i))));
}

double Gradient::relative_dist_energy(double dist, double l, double u) {
    if (dist > u) return square(dist - u); else if (dist < l) return square(l - dist); else return 0;
}

double Gradient::atom_dih_energy(int n) {
    double en = 0;
    for (auto && i : _involved_dihs[n]) {
        en += dih_row_energy(i);
    }
    return en;
}

double Gradient::dih_row_energy(const DihBound::key_type &ls) {
    double dih = geom::dihedral(c.row(ls[0]), c.row(ls[1]), c.row(ls[2]), c.row(ls[3]));
    return relative_dih_energy(dih, _dih_bound[ls]);
}

void Gradient::shrink_dihedral(double &n) {
    while (true) {
        if (n >= PI) n -= 2 * PI;
        else if (n < -PI) n += 2 * PI;
        else if (n < 0) n = -n;
        else break;
    }
}

double Gradient::relative_dih_energy(double dih, double std) {
    double d = dih - std; shrink_dihedral(d);
    return 10 * d * d;
}

} // namespace dg
END_JN

