#pragma once

#include "../pdb/Chain.hpp"
#include "../matrix.hpp"

namespace jian {

Chain psb2aa(const Eigen::MatrixXd &c, int beg, int end);

void psb_extract_frags(const std::string &s);

} // namespace jian

