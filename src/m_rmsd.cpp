#include <list>
#include <array>
#include <iostream>
#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/geom.hpp>
#include <jian/utils/string.hpp>

namespace jian {

	namespace {

		class RmsdComponent {
		public:
			using nums_t = std::vector<int>;

			bool m_loose;
			bool m_has_nums;
			bool m_has_traj;
			nums_t nums;
			Par::pars_t m_files;
			str_t m_ref;
			str_t m_tgt;

			RmsdComponent(const Par &par) {
				m_loose = par.has("loose");
				m_has_nums = par.has("nums");
				m_has_traj = par.has("traj");

				if (par.has("nums")) {
					set_nums(nums, par["nums"][0]);
				}

				par.setv(m_files, "s");
				par.set(m_ref, "ref");
				par.set(m_tgt, "traj");

			}

			void set_nums(nums_t &nums, const str_t &par) {
				std::vector<str_t> v, w;
				jian::tokenize(par, v, "+");
				for (auto && i : v) {
					jian::tokenize(i, w, "-");
					if (w.size() == 1) {
						nums.push_back(std::stoi(w[0]) - 1);
					}
					else if (w.size() == 2) {
						for (int j = std::stoi(w[0]); j <= std::stoi(w[1]); j++) {
							nums.push_back(j - 1);
						}
					}
				}
			}

			val_t rmsd(const Model &m1, const Model &m2) {
				assert(num_residues(m1) == num_residues(m2));
				int len = num_residues(m1);
				std::list<std::array<double, 3>> l1, l2;
				int i1 = 0, j1 = 0, i2 = 0, j2 = 0, n = 0;

				while (true) {
					const Residue & r1 = m1[i1][j1];
					const Residue & r2 = m2[i2][j2];
					if (!m_has_nums || std::find(nums.begin(), nums.end(), i1) != nums.end()) {
						if (!m_loose) assert(r1.name == r2.name);
						for (auto && a1 : r1) {
							for (auto && a2 : r2) {
								if (a1.name == a2.name) {
									l1.push_back({ a1[0],a1[1],a1[2] });
									l2.push_back({ a2[0],a2[1],a2[2] });
									n++;
								}
							}
						}
					}
					j1++;
					if (j1 >= m1[i1].size()) {
						j1 = 0;
						i1++;
						if (i1 >= m1.size()) {
							break;
						}
					}
					j2++;
					if (j2 >= m2[i2].size()) {
						j2 = 0;
						i2++;
						if (i2 >= m2.size()) {
							break;
						}
					}
				}

				Mat mat1(n, 3), mat2(n, 3);
				int i = 0;
				auto it1 = l1.begin(), it2 = l2.begin();
				for (; it1 != l1.end() && it2 != l2.end(); it1++, it2++, i++) {
					for (int k = 0; k < 3; k++) {
						mat1(i, k) = (*it1)[k];
						mat2(i, k) = (*it2)[k];
					}
				}
				return geom::rmsd(mat1, mat2);
			}

			void run() {
				if (m_has_traj) {
					Model ref;
					mol_read(ref, m_ref);
					for_each_model(m_tgt, [this, &ref](const Model &m, int i) {
						JN_OUT << m.num << ' ' << rmsd(ref, m) << std::endl;
					});
				}
				else {
					JN_OUT << rmsd(mol_read_to<Model>(m_files[0]), mol_read_to<Model>(m_files[1])) << std::endl;
				}
			}
		};

		REGISTER_NSP_COMPONENT(rmsd) {
			RmsdComponent rmsd(par);
			rmsd.run();
		}
	}

} // namespace jian

