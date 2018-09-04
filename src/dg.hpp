#pragma once

#include "jian.hpp"

namespace jian {

struct DgConstraint {
    std::vector<int> key;
    std::vector<double> val;
};

class DgConstraints : public std::vector<DgConstraint> {
public:
    template<typename K1, typename K2>
    bool key_equal(const K1 &k1, const K2 &k2) {
        if (k1.size() != k2.size()) return false;
        return std::all_of(k1.begin(), k1.end(), [&k2](auto &&n1){
            return std::any_of(k2.begin(), k2.end(), [&n1](auto &&n2){
                return n1 == n2;
            });
        });
    }

    void append(const DgConstraint &constraint) {
        for (auto && c : *this) {
            if (key_equal(c.key, constraint.key)) {
                c.val = constraint.val;
                return;
            }
        }
        this->push_back(constraint);
    }
};

class DgImpl;

class DG {
    DgImpl *impl;
public:
    DG();

    void read(int len_,
              const std::array<double, 2> &default_range_,
              const std::vector<DgConstraint> &distance_constraints_,
              const std::vector<DgConstraint> &dihedral_constraints_ = {});

    Mat sample();

    void log_file(std::string);

    ~DG();

};

}

