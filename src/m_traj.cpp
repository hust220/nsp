#include <list>
#include <numeric>
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
			str_t m_func;
			str_t m_traj;
			str_t m_ref_file;
			str_t m_align;
			Chain m_ref;
			std::deque<int> m_common_ref, m_common_tgt;
			int m_bin = 1;
			bool m_loose;
			Mat m_mat_ref;
			Mat m_mat_tgt;

			TrajComponent(const Par &par) {
				auto v = par.getv("global");
				assert(v.size() >= 3);
				m_func = v[1];
				to_lower(m_func);

				m_traj = v[2];
				if (v.size() == 4) {
					m_ref_file = v[3];
					chain_read_model(m_ref, m_ref_file);
				}

				par.set(m_bin, "b", "bin");

				m_loose = par.has("loose");

				par.set(m_align, "align");
				set_common();

				//set_mat(m_mat_ref, m_ref, m_common_ref);
			}

			//void set_mat(Mat &mat, const Chain &chain, const std::deque<int> &common) {
			//	int i, j, n, l, rows;

			//	l = size(chain);
			//	rows = size(common);
			//	mat.resize(rows, 3);
			//	n = 0;
			//	for (i = 0; i < l; i++) {
			//		if (std::find(common.begin(), common.end(), i) != common.end()) {
			//			const Atom &atom = chain[i]["C4*"];
			//			for (j = 0; j < 3; j++) {
			//				mat(n, j) = atom[j];
			//			}
			//			n++;
			//		}
			//	}
			//}

			Chain model_to_chain(const Model &m) {
				Chain c;
				for (auto && chain : m) {
					for (auto && res : chain) {
						c.push_back(res);
					}
				}
				return c;
			}

			val_t rmsd(const Chain &tgt) {
				std::deque<std::array<double, 3>> l1, l2;
				int i, l, n;

				auto it_ref = m_common_ref.begin();
				auto it_tgt = m_common_tgt.begin();
				while (it_ref != m_common_ref.end() && it_tgt != m_common_tgt.end()) {
					const Residue & r1 = m_ref[*it_ref];
					const Residue & r2 = tgt[*it_tgt];
					if (!m_loose) assert(r1.name == r2.name);
					for (auto && a1 : r1) {
						for (auto && a2 : r2) {
							if (a1.name == a2.name) {
								l1.push_back({ a1[0],a1[1],a1[2] });
								l2.push_back({ a2[0],a2[1],a2[2] });
							}
						}
					}
					it_ref++;
					it_tgt++;
				}
				

				i = 0;
				l = l1.size();
				Mat mat1(l, 3);
				Mat mat2(l, 3);
				for (i = 0; i < l; i++) {
					for (int k = 0; k < 3; k++) {
						mat1(i, k) = l1[i][k];
						mat2(i, k) = l2[i][k];
					}
				}
				return geom::rmsd(mat1, mat2);
			}

			void set_common() {
				std::ifstream ifile;
				str_t seq1, seq2;
				std::deque<int> dq1, dq2;
				int i, l, n1, n2;

				if (m_align.empty()) {
					l = size(m_ref);
					m_common_ref.resize(l);
					m_common_tgt.resize(l);
					std::iota(m_common_ref.begin(), m_common_ref.end(), 0);
					std::iota(m_common_tgt.begin(), m_common_tgt.end(), 0);
				}
				else {
					FOPEN(ifile, m_align);
					ifile >> seq1 >> seq2;
					assert(size(seq1) == size(seq2));
					l = size(seq1);
					n1 = 0;
					n2 = 0;
					for (i = 0; i < l; i++) {
						if (seq1[i] != '-' && seq2[i] != '-') {
							dq1.push_back(n1);
							dq2.push_back(n2);
						}
						if (seq1[i] != '-') n1++;
						if (seq2[i] != '-') n2++;
					}
					FCLOSE(ifile);
					m_common_ref = std::move(dq1);
					m_common_tgt = std::move(dq2);
				}
			}

			void rmsd() {
				for_each_model(m_traj, [this](const Model &model, int i) {
					if (i % m_bin == 0) {
						JN_OUT << i+1 << ' ' << rmsd(model_to_chain(model)) << std::endl;
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

