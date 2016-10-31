#include <vector>
#include <string>
#include <set>
#include <array>
#include <algorithm>
#include "CG6p.hpp"

namespace jian {
	REG_CG("6p", CG6p);

	CG6p::CG6p() {
		m_cg = "6p";
	}

	int CG6p::res_size() const {
		return 6;
	}

	const std::vector<std::string> CG6p::m_basic_atoms{
		"C5*", "C1*", "O3*", "C2", "C4", "C6"
	};

	Residue CG6p::to_cg(const Residue &r) const {
		if (is_cg(r)) {
			return r;
		}
		else {
			Residue res;
			res.name = r.name;
			res.m_cg = m_cg;
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

} // namespace jian
