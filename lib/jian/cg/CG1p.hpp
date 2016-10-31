#pragma once

#include "../matrix.hpp"
#include "../pdb/Chain.hpp"
#include "CG.hpp"

namespace jian {

	class CG1p : public CG {
	public:
		CG1p();
		virtual Residue to_cg(const Residue &r) const;
		virtual int res_size() const;
	};

} // namespace jian

