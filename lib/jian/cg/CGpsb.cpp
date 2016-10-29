#include <vector>
#include <string>
#include <set>
#include <array>
#include <algorithm>
#include "CGpsb.hpp"

namespace jian {

	Residue CGpsb::res(const Residue &r) {
		static std::vector<std::string> v{ "N1", "C2", "N3", "C4", "C5", "C6" };
		if (r.is_cg<CGpsb>()) {
			return r;
		}
		else {
			Residue res;
			res.name = r.name;
			std::array<double, 3> arr{ 0, 0, 0 };
			Atom atm;
			for (auto && atom : r) {
				if (atom.name == "C4*") {
					atm.init("P", atom[0], atom[1], atom[2]);
					res.push_back(atm);
				}
				else if (atom.name == "C1*") {
					atm.init("S", atom[0], atom[1], atom[2]);
					res.push_back(atm);
				}
				else if (std::find(v.begin(), v.end(), atom.name) != v.end()) {
					for (int i = 0; i < 3; i++) {
						arr[i] += atom[i];
					}
				}
			}
			for (int i = 0; i < 3; i++) {
				arr[i] /= 6.0;
			}
			atm.init("B", arr[0], arr[1], arr[2]);
			res.push_back(atm);
			return res;
		}
	}

	//Chain CGpsb::chain(const Chain &chain) {
	//    Chain c;
	//    c.name = chain.name;
	//    c.model_name = chain.model_name;
	//    for (auto && r : chain) {
	//		Residue res = r;
	//        c.push_back(res.cg<CGpsb>());
	//    }
	//    return c;
	//}

} // namespace jian
