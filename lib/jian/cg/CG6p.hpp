#pragma once

#include "../matrix.hpp"
#include "../pdb/Chain.hpp"
#include "CG.hpp"

namespace jian {

	class CG6p : public CG<CG6p> {
	public:
		static const int size_res;
		static const Residue::cg_code m_cg;
		static const std::vector<std::string> m_basic_atoms;

		static Residue res(const Residue &r);
		//static Chain chain(const Chain &chain);
		static Chain aa(const Mat &c, int beg, int end);
		static void extract_frags(const std::string &s);
	};

} // namespace jian


