#pragma once

#include "../matrix.hpp"
#include "../pdb/Chain.hpp"

namespace jian {

class CGpsb {
public:
    static const int size_res = 3;
	static const Residue::cg_code m_cg = Residue::CG_PSB;

    static Residue res(const Residue &r);
    static Chain chain(const Chain &chain);
    static Chain aa(const Mat &c, int beg, int end);
    static void extract_frags(const std::string &s);
    static bool is_psb(const Residue &r);
};

} // namespace jian


