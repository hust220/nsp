#include <vector>
#include <string>
#include <set>
#include <array>
#include <algorithm>
#include "CG6p.hpp"

namespace jian {
	const int CG6p::size_res = 5;

	const Residue::cg_code CG6p::m_cg = Residue::CG_6P;

	const std::vector<std::string> CG6p::m_basic_atoms{
		"C5*", "C1*", "O3*", "C2", "C4", "C6"
	};

	Residue CG6p::res(const Residue &r) {
		if (r.is_cg<CG6p>()) {
			return r;
		}
		else {
			Residue res;
			res.name = r.name;
			for (auto && name : m_basic_atoms) {
				auto it = std::find_if(r.begin(), r.end(), [&name](const Atom &atom) {
					return atom.name == name;
				});
				if (it == r.end()) {
					throw "CG6p::res atom '" + name + "' not found!";
				}
				else {
					res.push_back(*it);
				}
			}
			return res;
		}
	}

	//Chain CG6p::chain(const Chain &chain) {
	//	Chain c;
	//	c.name = chain.name;
	//	c.model_name = chain.model_name;
	//	for (auto && r : chain) {
	//		Residue res = r;
	//		c.push_back(res.cg<CG6p>());
	//	}
	//	return c;
	//}

} // namespace jian
