#include "../geom.hpp"
#include "../utils/Env.hpp"
#include "CompleteResidue.hpp"

namespace jian {

	CompleteResidue::CompleteResidue() {
		std::string path = Env::lib() + "/RNA/pars/cg/CompleteResidue/";
		for (std::string s : {"RNA", "DNA"}) {
			mol_read(m_bb[s], path + s + ".bb.pdb");
			const pdb::Names &names = pdb::Names::instance(s);
			for (auto && name : names.res) {
				mol_read(m_base[name], path + name + ".base.pdb");
			}
		}
	}

	const CompleteResidue &CompleteResidue::instance() {
		static const CompleteResidue complete;
		return complete;
	}

	bool CompleteResidue::lack_atoms(const Residue &res) const {
		const pdb::names_t &names = pdb::res_included_atoms(res.name);
		return std::any_of(names.begin(), names.end(), [&res](std::string name) {
			return std::none_of(res.begin(), res.end(), [&name](const Atom &atom) {
				return atom.name == name;
			});
		});
	}

	void CompleteResidue::operator()(Residue &res) const {
		auto align = [](Residue &r1, const Residue &r2) {
			int l = 0;
			std::deque<std::array<double, 3>> dq1, dq2;
			for (auto && atom1 : r1) {
				auto it = std::find_if(r2.begin(), r2.end(), [&atom1](const Atom &atom2) {
					return atom1.name == atom2.name;
				});
				if (it != r2.end()) {
					dq1.push_back({ atom1[0], atom1[1], atom1[2] });
					dq2.push_back({ it->at(0), it->at(1), it->at(2) });
					l++;
				}
			}

			Mat m1(l, 3), m2(l, 3);
			for (int i = 0; i < l; i++) {
				for (int j = 0; j < 3; j++) {
					m1(i, j) = dq1[i][j];
					m2(i, j) = dq2[i][j];
				}
			}

			geom::Superposition<double> sp(m1, m2);
			for (auto && atom : r1) {
				sp.apply(atom);
			}
		};

		if (lack_atoms(res)) {
			Residue bb = m_bb.at(pdb::res_mol_type<std::string>(res.name));
			Residue base = m_base.at(res.name);

			align(bb, res);
			align(base, res);

			res.clear();
			res += bb;
			res += base;
		}

	}
}