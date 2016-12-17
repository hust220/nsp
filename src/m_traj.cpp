#include <list>
#include <numeric>
#include <array>
#include <iostream>
#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/geom.hpp>
#include <jian/utils/string.hpp>

BEGIN_JN

	namespace {

		class TrajComponent {
		public:
			Str m_func;
			Str m_traj;
			Str m_ref_file;
			Str m_align;
			Chain m_ref;
			std::deque<int> m_common_ref, m_common_tgt;
			int m_bin = 1;
			bool m_loose;
			Mat m_mat_ref;
			Mat m_mat_tgt;
			Par m_par;

			TrajComponent(const Par &par) {
				m_par = par;

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

			}

			template<typename _Chain>
			Num rmsd(const _Chain &tgt) {
				std::deque<std::array<double, 3>> l1, l2;
				int i, l;

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
				Str seq1, seq2;
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
						JN_OUT << i+1 << ' ' << rmsd(model.residues()) << std::endl;
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

			void extract() {
				Num rmsd_target = -1;
				Num num_target = -1;
				m_par.set(rmsd_target, "rmsd");
				m_par.set(num_target, "n");
				MolReader mol_reader(m_traj);
				auto it = mol_reader.model_begin();
				auto it_target = it;
				Num min = 999;
				for (; it != mol_reader.model_end(); it++)
				{
					if (num_target != -1)
					{
						if (it.n + 1 == num_target)
						{
							it_target = it;
							JN_OUT << *it;
							break;
						}
					}
					else
					{
						Num d = STD_ fabs(rmsd_target - rmsd(it->residues()));
						if (d < min)
						{
							it_target = it;
							min = d;
						}
					}
				}
				JN_OUT << *it_target;
			}

			void run() {
				if (m_func == "rmsd") {
					rmsd();
				}
				else if (m_func == "compress") {
					compress();
				}
				else if (m_func == "extract") {
					extract();
				}
			}
		};

		REGISTER_NSP_COMPONENT(traj) {
			TrajComponent traj(par);
			traj.run();
		}
	}

END_JN

