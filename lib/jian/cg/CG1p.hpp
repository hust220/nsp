#pragma once

#include "../matrix.hpp"
#include "../pdb/Chain.hpp"

namespace jian {

struct CG1p {
    static int size_res;
    static Residue res(const Residue &r);
    static Chain chain(const Chain &chain);
    static Chain aa(const Mat &c, int beg, int end);
};

} // namespace jian

