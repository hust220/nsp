#pragma once

#include "../matrix.hpp"
#include "../pdb/Chain.hpp"
#include "CG.hpp"

BEGIN_JN

	class CG6p : public CG {
	public:
		static const pdb::names_t m_basic_atoms;

		CG6p();
		//Residue get_atoms(const Residue &res, const pdb::names_t &names) const;
		virtual Residue to_cg(const Residue &r) const;
		virtual int res_size() const;
	};

END_JN


