#pragma once

#include "../matrix.hpp"
#include "../pdb/Chain.hpp"
#include "CG.hpp"

BEGIN_JN

	class CGaa : public CG {
	public:
		static const std::vector<std::string> m_basic_atoms;

		CGaa();
		virtual Residue to_cg(const Residue &r) const;
		virtual int res_size() const;
	};

END_JN


