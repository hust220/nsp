#pragma once

#include "../pdb.hpp"
#include "../matrix.hpp"

namespace jian {

struct scoring {
    static double new_score(const Model &model);
    static double new_score(const Eigen::MatrixXd &mat, int index);
};

}

