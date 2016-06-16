#pragma once

#include "../pdb/Chain.hpp"
#include "../matrix.hpp"

namespace jian {

Chain c2a(const Eigen::MatrixXd &c, int beg, int end);

} // namespace jian

