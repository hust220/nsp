#include <vector>
#include <string>
#include <set>
#include <array>
#include <algorithm>
#include "cg_psb.hpp"

BEGIN_JN
	REG_CG("psb", CGpsb);

	CGpsb::CGpsb() {
		m_cg = "psb";
	}

	int CGpsb::res_size() const {
		return 3;
	}

	Residue CGpsb::to_cg(const Residue &r) const {
		static std::vector<std::string> v{ "N1", "C2", "N3", "C4", "C5", "C6" };
		if (is_cg(r)) {
			return r;
		}
		else {
			Residue res;
			res.name = r.name;
			res.m_cg = m_cg;
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

END_JN
