#include "CG1p.hpp"

namespace jian {

	REG_CG("1p", CG1p);

	CG1p::CG1p() {
		m_cg = "1p";
	}

	int CG1p::res_size() const {
		return 1;
	}

	Residue CG1p::to_cg(const Residue &r) const {
		Residue res;
		res.name = r.name;
		res.m_cg = m_cg;
		for (auto && atom : r) {
			if (atom.name == "C4*") {
				res.push_back(atom);
			}
		}
		return res;
	}

} // namespace jian
