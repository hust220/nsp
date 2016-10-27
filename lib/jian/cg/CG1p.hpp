#pragma once

#include "../matrix.hpp"
#include "../pdb/Chain.hpp"

namespace jian {

struct CG1p {
	static const Residue::cg_code m_cg = Residue::CG_1P;
    static const int size_res = 1;

    static Residue res(const Residue &r);
    static Chain chain(const Chain &chain);
    static Chain aa(const Mat &c, int beg, int end);
};

} // namespace jian

