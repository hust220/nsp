#include <list>
#include <array>
#include <iostream>
#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/geom.hpp>
#include <jian/utils/string.hpp>

namespace jian {

	namespace {

		class TrajComponent {
		public:
			std::string m_func;
			std::string m_traj;
			std::string m_ref;
			int m_bin = 1;
			bool m_loose;

			TrajComponent(const Par &par) {
				auto v = par.getv("global");
				assert(v.size() >= 3);
				m_func = v[1];
				to_lower(m_func);

				m_traj = v[2];
				if (v.size() == 4) {
					m_ref = v[3];
				}
				par.set(m_bin, "b", "bin");
				m_loose = par.has("loose");
			}

			val_t rmsd(const Model &m1, const Model &m2) {
				assert(num_residues(m1) == num_residues(m2));
				int len = num_residues(m1);
				std::list<std::array<double, 3>> l1, l2;
				int i1 = 0, j1 = 0, i2 = 0, j2 = 0, n = 0;

				while (true) {
					const Residue & r1 = m1[i1][j1];
					const Residue & r2 = m2[i2][j2];
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

			void rmsd() {
				Model m;
				mol_read(m, m_ref);
				for_each_model(m_traj, [this, &m](const Model &model, int i) {
					if (i % m_bin == 0) {
						JN_OUT << i+1 << ' ' << rmsd(m, model) << std::endl;
					}
				});
			}

			void compress() {
				MolWriter writer(JN_OUT);
				for_each_model(m_traj, [this, &writer](const Model &model, int i) {
					if (i % m_bin == 0) {
						LOG << "Reading: model-" << i + 1 << std::endl;
						writer.write(model);
					}
				});
			}

			void run() {
				if (m_func == "rmsd") {
					rmsd();
				}
				else if (m_func == "compress") {
					compress();
				}
			}
		};

		REGISTER_NSP_COMPONENT(traj) {
			TrajComponent traj(par);
			traj.run();
		}
	}

} // namespace jian

