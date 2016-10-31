#pragma once

#include "../matrix.hpp"
#include "../pdb/Chain.hpp"
#include "CG.hpp"

namespace jian {

	class CG6p : public CG {
	public:
		static const std::vector<std::string> m_basic_atoms;

		CG6p();
		virtual Residue to_cg(const Residue &r) const;
		virtual int res_size() const;
	};

} // namespace jian


