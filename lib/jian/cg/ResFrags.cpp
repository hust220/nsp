#include "ResFrags.hpp"
#include "CG.hpp"
#include "../utils/Env.hpp"
#include "../geom.hpp"

BEGIN_JN

std::map<std::string, std::map<int, ResFrags>> ResFrags::m_instances;

ResFrags &ResFrags::instance(S cg, int frag_size) {
	if (m_instances.find(cg) == m_instances.end() || m_instances[cg].find(frag_size) == m_instances[cg].end()) {
		m_instances[cg][frag_size] = ResFrags();
		m_instances[cg][frag_size].init(cg, frag_size);
	}
	return m_instances[cg][frag_size];
}

ResFrags::~ResFrags() {
	for (auto && frag : m_mats) {
		delete frag;
	}
	for (auto && chain : m_chains_aa) {
		delete chain;
	}
	for (auto && chain : m_chains_cg) {
		delete chain;
	}
	delete m_cg;
}

void ResFrags::init(S cg, int frag_size) {
	m_cg = CG::fac_t::create(cg);
	m_res_size = m_cg->res_size();
	m_frag_size = frag_size;
	m_path = Env::lib() + "/RNA/pars/cg/CG2AA/templates.pdb";
	extract_frags(m_path);
}

Chain ResFrags::get_chain(int i, const Mat &c) {
	geom::Superposition<double> sp;
	Chain chain;

	//chain = m_chains[i];
	//sp.init(*(m_frags[i]), c);
	//for (auto && res : chain) for (auto && atom : res) {
	//	sp.apply(atom);
	//}
	return chain;
}

void ResFrags::extract_frags(const S &pdb) {
	std::deque<int> dq;
	Chain full_chain, chain;
	Chain *p_chain_aa, *p_chain_cg;
	Mat *p_mat;
	int i, j, k, n;
	S name;

	chain_read_model(full_chain, pdb);
	chain = m_cg->to_cg(full_chain);
	for (n = 0; n < chain.size(); n++) {
		if (dq.size() >= 1 && geom::distance(chain[dq.back()][2], chain[n][0]) > 4) {
			dq.clear();
		}
		dq.push_back(n);
		if (dq.size() == m_frag_size) {
			std::ostringstream stream;
			p_mat = new Mat(m_frag_size * m_res_size, 3);
			m_mats.push_back(p_mat);
			p_chain_aa = new Chain;
			m_chains_aa.push_back(p_chain_aa);
			p_chain_cg = new Chain;
			m_chains_cg.push_back(p_chain_cg);
			for (i = 0; i < m_frag_size; i++) {
				stream << full_chain[dq[i]].name;
				p_chain_aa->push_back(full_chain[dq[i]]);
				p_chain_cg->push_back(chain[dq[i]]);
				for (j = 0; j < m_res_size; j++) {
					for (k = 0; k < 3; k++) {
						(*p_mat)(i*m_res_size + j, k) = chain[dq[i]][j][k];
					}
				}
			}
			//*p_chain_cg = m_cg->to_cg(*p_chain_aa);
			m_ids[stream.str()].push_back(m_chains_aa.size() - 1);
			dq.pop_front();
		}
	}
}

END_JN