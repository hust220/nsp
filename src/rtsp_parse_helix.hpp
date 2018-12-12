#pragma once

#include <Eigen/Dense>
#include "pdb.hpp"

namespace jian {

struct parse_helix_t {
    using Vec = Eigen::Matrix<double, 3, 1>;
    Vec origin, x, y, z; 
    double theta, phi;
};

parse_helix_t parse_helix(const Model &helix);
parse_helix_t parse_helix(const Mat &helix);
Eigen::MatrixXd make_standard_helix(int);
void dihs_std_helix();

}

