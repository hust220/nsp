#pragma once

#include "jian/matrix.hpp"
#include "../pdb/Chain.hpp"
#include "CG.hpp"

BEGIN_JN

	class CG1p : public CG {
	public:
		CG1p();
		virtual Residue to_cg(const Residue &r) const;
		virtual int res_size() const;
	};

END_JN

