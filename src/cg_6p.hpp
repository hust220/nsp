#pragma once

#include "matrix.hpp"
#include "pdb_chain.hpp"
#include "cg.hpp"

namespace jian {

	class CG6p : public CG {
	public:
		static const pdb::names_t m_basic_atoms;

		CG6p();
		//Residue get_atoms(const Residue &res, const pdb::names_t &names) const;
		virtual Residue to_cg(const Residue &r) const;
		virtual int res_size() const;
	};

}


