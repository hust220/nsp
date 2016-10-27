#include <string>
#include "../utils/Env.hpp"
#include "../utils/Par.hpp"
#include "../pdb.hpp"
#include "ScoreAa.hpp"

namespace jian {

	REG_SCORER("aa", Score<AA>);

	void Score<AA>::init() {
		std::string lib = Env::lib() + "/RNA/pars/scoring";

		m_bin_dist = 0.3;
		m_bin_dih = 4.5;
		m_cutoff = 20;

		m_score_dist = 0;
		m_score_dih = 0;

		m_constant = 27.1118;
		m_weight_dist = 0.433513;
		m_weight_dih = 1.59348;

		m_par_dist = lib + "/dist.frequencies";
		m_par_dih = lib + "/par_dih";

		m_dist_anal = new DistAnal;
		m_dist_anal->init(m_bin_dist, m_cutoff);
		m_dist_anal->read_freqs(m_par_dist);

		m_dih_anal = new DihAnal(m_bin_dih);
		m_dih_anal->read_parm(m_par_dih);
	}

	Score<AA>::~Score() {
		if (m_dist_anal != NULL) delete m_dist_anal;
		if (m_dih_anal != NULL) delete m_dih_anal;
	}

	void Score<AA>::run(const Chain &chain) {
		m_score_dist = m_dist_anal->run(chain).score;
		m_score_dih = m_dih_anal->run(chain).score;
		m_score = m_constant + m_score_dist * m_weight_dist + m_score_dih * m_weight_dih;
	}

	void Score<AA>::train(const Chain &c) {
		m_dist_anal->train(c);
	}

	void Score<AA>::print_counts(std::ostream & stream) const {
		m_dist_anal->print_counts(stream);
	}

	void Score<AA>::print_freqs(std::ostream & stream) const {
		m_dist_anal->print_freqs(stream);
	}

	double Score<AA>::en_stacking(const Residue &r1, const Residue &r2) {
		return m_dist_anal->en_stacking(r1, r2);
	}

	double Score<AA>::en_pairing(const Residue &r1, const Residue &r2) {
		return m_dist_anal->en_pairing(r1, r2);

	}


} // namespace jian

