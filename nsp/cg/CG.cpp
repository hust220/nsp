#include <deque>
#include "jian/matrix.hpp"
#include "../pdb.hpp"
#include "CG.hpp"
#include "jian/utils/Env.hpp"
#include "jian/utils/log.hpp"
#include "jian/geom.hpp"
#include "../nuc3d/Convert.hpp"

BEGIN_JN
	Chain CG::to_cg(const Chain &chain) const {
		Chain c;
		c.name = chain.name;
		c.model_name = chain.model_name;
		c.m_cg = m_cg;
		for (auto && r : chain) {
			c.push_back(to_cg(r));
		}
		return c;
	}

	Model CG::to_cg(const Model &model) const {
		Model m;
		m.name = model.name;
		m.type = model.type;
		m.num = model.num;
		m.m_cg = m_cg;
		for (auto && c : model) {
			m.push_back(to_cg(c));
		}
		return m;
	}

	Molecule CG::to_cg(const Molecule &molecule) const {
		Molecule mol;
		mol.name = molecule.name;
		mol.m_cg = m_cg;
		for (auto && m : molecule) {
			mol.push_back(to_cg(m));
		}
		return mol;
	}

	class CG2AA {
	public:
		CG *m_cg;
		std::deque<Mat *> m_frags;
		std::deque<Chain> m_chains;
		std::deque<std::string> m_names;
		S m_path;
		int m_frag_size = 4;
		int m_res_size;
		static std::map<std::string, CG2AA> m_instances;

		static CG2AA &instance(S cg) {
			if (m_instances.find(cg) == m_instances.end()) {
				m_instances[cg] = CG2AA();
				m_instances[cg].init(cg);
			}
			return m_instances[cg];
		}

		~CG2AA() {
			for (auto && frag : m_frags) {
				delete frag;
			}
			delete m_cg;
		}

		void init(S cg) {
			m_cg = CG::fac_t::create(cg);
			m_res_size = m_cg->res_size();
			m_path = Env::lib() + "/RNA/pars/cg/CG2AA/templates.pdb";
			extract_frags(m_path);
		}

		Chain get_chain(int i, const Mat &c) {
			geom::Superposition<double> sp;
			Chain chain;
			
			chain = m_chains[i];
			sp.init(*(m_frags[i]), c);
			for (auto && res : chain) for (auto && atom : res) {
				sp.apply(atom);
			}
			return chain;
		}

		template<typename Coord, typename Frag>
		auto run(Coord &&coord, Frag &&frag) {
			int i, j, k, num_atoms, len, index;
			Chain residues;

			//LOG << "## CG2AA\n" << std::endl;
			num_atoms = m_frag_size * m_res_size;
			len = frag[1] - frag[0] + 1;
			Chain chain;
			for (i = 0; i < len - num_atoms + 1; i += m_res_size) {
				Mat c(num_atoms, 3);
				for (j = 0; j < num_atoms; j++) {
					for (k = 0; k < 3; k++) {
						c(j, k) = coord(frag[0] + i + j, k);
					}
				}
				std::deque<double> scores;
				for (auto && frag : m_frags) {
					scores.push_back(geom::rmsd(*frag, c));
				}
				auto min = std::min_element(scores.begin(), scores.end());
				index = std::distance(scores.begin(), min);
				residues = get_chain(index, c);
				if (i == 0) {
					for (j = 0; j < m_frag_size - 1; j++) {
						chain.push_back(residues[j]);
					}
				}
				chain.push_back(residues[m_frag_size - 1]);
			}
			return chain;
		}

		void extract_frags(const S &pdb) {
			std::deque<int> dq;
			Chain full_chain, chain;
			int i, j, k, n;
			S name;

			chain_read_model(full_chain, pdb);
			chain = m_cg->to_cg(full_chain);
			for (n = 0; n < chain.size(); n++) {
				if (dq.size() >= 1 && geom::distance(chain[dq.back()][0], chain[n][0]) > 10) {
					dq.clear();
				}
				dq.push_back(n);
				if (dq.size() == m_frag_size) {
					m_chains.push_back(Chain{});
					Mat *m = new Mat(m_frag_size * m_res_size, 3);
					for (i = 0; i < m_frag_size; i++) {
						m_chains.back().push_back(full_chain[dq[i]]);
						for (j = 0; j < m_res_size; j++) {
							for (k = 0; k < 3; k++) {
								(*m)(i*m_res_size + j, k) = chain[dq[i]][j][k];
							}
						}
					}
					m_frags.push_back(m);
					dq.pop_front();
				}
			}
		}

	};

	std::map<std::string, CG2AA> CG2AA::m_instances;

	Mat CG::chain_to_coords(const Chain &chain) const {
		int n_atoms = num_atoms(chain);
		Mat c(n_atoms, 3);
		int n_atom = 0;
		for (int i = 0; i < size(chain); i++) {
			for (auto && atom : chain[i]) {
				mat_set_rows(c, n_atom, atom);
				n_atom++;
			}
		}
		return c;
	}

	Chain CG::to_aa(const Chain &chain) const {
		Mat &&mat = chain_to_coords(chain);
		Str s = seq(chain);
		Chain &&c = to_aa(mat, 0, mat.rows() - 1);
		int i = 0;
		for (Residue &res : c) {
			res = convert_res(res, to_str(s[i]));
			i++;
		}
		return c;
	}

    Model CG::to_aa(const Model &model) const {
        Model m;
        m.name = model.name;
        m.type = model.type;
        m.num = model.num;
        for (auto && chain : model) {
            m.push_back(to_aa(chain));
        }
        m.m_cg = "aa";
        return m;
    }

	Chain CG::to_aa(const Mat &c, int beg, int end) const {
		return CG2AA::instance(m_cg).run(c, std::vector<int>{beg, end});
	}


END_JN
