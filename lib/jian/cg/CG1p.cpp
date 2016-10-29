#include "CG1p.hpp"

namespace jian {

	//Residue::cg_code CG1p::m_cg = Residue::CG_1P;

	Residue CG1p::res(const Residue &r) {
		Residue res;
		res.name = r.name;
		for (auto && atom : r) {
			if (atom.name == "C4*") {
				res.push_back(atom);
			}
		}
		return res;
	}

	//Chain CG1p::chain(const Chain &chain) {
	//    Chain c;
	//    c.name = chain.name;
	//    c.model_name = chain.model_name;
	//    for (auto && r : chain) {
	//		Residue res = r;
	//        c.push_back(res.cg<CG1p>());
	//    }
	//    return c;
	//}

} // namespace jian
