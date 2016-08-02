#pragma once

#include "../pdb.hpp"
#include "../matrix.hpp"

namespace jian {
namespace scoring {

double new_score(const Model &model);
double new_score(const Eigen::MatrixXd &mat, int index);

}
}

