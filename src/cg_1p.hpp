#pragma once

#include "matrix.hpp"
#include "pdb_chain.hpp"
#include "cg.hpp"

BEGIN_JN

	class CG1p : public CG {
	public:
		CG1p();
		virtual Residue to_cg(const Residue &r) const;
		virtual int res_size() const;
	};

END_JN

