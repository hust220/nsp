#pragma once

#include "Job.hpp"
#include <vector>
#include <list>

namespace jian {
namespace dg {

class Gradient : public Job {
public:
    std::vector<std::list<DihEntry>> _involved_dihs;

    Gradient();
    Gradient(const Gradient &) = default;
    Gradient &operator =(const Gradient &) = default;

    void gradient();
    void init_dihs();
    double total_energy(const Mat &coord);
    double total_dist_energy(const Mat &coord);
    double total_dih_energy();
    double atom_energy(const Mat &coord, int n);
    double atom_dist_energy(const Mat &coord, int n);
    double atom_pair_energy(const Mat &coord, int i, int j);
    double relative_dist_energy(double dist, double l, double u);
    double atom_dih_energy(int n);
    double dih_row_energy(const DihBound::key_type &ls);
    void shrink_dihedral(double &n);
    double relative_dih_energy(double dih, double std);
};

} // namespace dg
} // namespace jian

