#pragma once

#include "../pdb/Chain.hpp"
#include "../matrix.hpp"

namespace jian {

Chain cg2aa(const Eigen::MatrixXd &c, int beg, int end);

} // namespace jian

