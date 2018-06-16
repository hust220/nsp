#include <vector>
#include <string>
#include <set>
#include <array>
#include <algorithm>
#include "cg_6p.hpp"
#include "cg_complete_res.hpp"

namespace jian {
	REG_CG("6p", CG6p);

	CG6p::CG6p() {
		m_cg = "6p";
	}

	int CG6p::res_size() const {
		return 6;
	}

	const std::vector<std::string> CG6p::m_basic_atoms{
		"P", "C4*", "C1*", "C2", "C4", "C6"
	};

	//Residue CG6p::get_atoms(const Residue &res, const pdb::names_t &names) const {
	//	Residue r;
	//	r.m_cg = res.m_cg;
	//	r.name = res.name;
	//	for (auto && atom : res) {
	//		if (std::find(names.begin(), names.end(), atom.name) != names.end()) {
	//			r.push_back(atom);
	//		}
	//	}
	//	return r;
	//}

	Residue CG6p::to_cg(const Residue &r) const {
		auto foo = [this](const Residue &res, auto && names) {
			Residue r;
			r.m_cg = m_cg;
			r.name = res.name;
			for (auto &&name : names) {
				auto it = std::find_if(res.begin(), res.end(), [&name](auto &&atom) {
					return atom.name == name;
				});
				std::ostringstream stream;
				stream << "jian::CG6p::to_cg error! Atom '" << name << "' not found";
				if (it == res.end()) throw stream.str();
				r.push_back(*it);
			}
			return r;
		};

		const CompleteResidue &complete = CompleteResidue::instance();

		if (is_cg(r)) {
			return r;
		}
		else {
			if (complete.lack_atoms(r)) {
				Residue res = r;
				complete(res);
				return foo(res, m_basic_atoms);
			}
			else {
				return foo(r, m_basic_atoms);
			}

		}
	}

}
