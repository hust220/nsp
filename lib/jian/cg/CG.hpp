#pragma once

#include <map>
#include <string>
#include "../pdb.hpp"
#include "../utils/Factory.hpp"

#define REG_CG(name, type) REGISTER_FACTORY(jian::CG::constructor_t, name, type)

BEGIN_JN

	class CG {
	public:
		using constructor_t = CG*(void);
		using fac_t = Factory<CG::constructor_t>;

		Str m_cg;
		//std::map<Str, std::map<Str, std::map<Str, std::vector<Str>>>>
		//	m_maps;

		template<typename T>
		bool is_cg(T &&t) const {
			return t.m_cg == m_cg;
		}

		virtual Residue to_cg(const Residue &residue) const = 0;
		virtual int res_size() const = 0;

		Chain to_cg(const Chain &chain) const;

		Model to_cg(const Model &model) const;

		Molecule to_cg(const Molecule &molecule) const;

		Mat chain_to_coords(const Chain &chain) const;

		Chain to_aa(const Chain &chain) const;

		Chain to_aa(const Mat &c, int beg, int end) const;

	};

}
