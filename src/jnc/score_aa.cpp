#include <string>
#include "env.hpp"
#include "par.hpp"
#include "pdb.hpp"
#include "score_aa.hpp"

namespace jian {

REG_SCORER("aa", ScoreAa);

ScoreAa::ScoreAa() {
	m_cg = CG::fac_t::create("aa");
}

void ScoreAa::init() {
	S lib = Env::lib() + "/RNA/pars/scoring";

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

ScoreAa::~ScoreAa() {
	if (m_dist_anal != NULL) delete m_dist_anal;
	if (m_dih_anal != NULL) delete m_dih_anal;
}

void ScoreAa::run(const Chain &chain) {
	m_score_dist = m_dist_anal->run(chain).score;
	m_score_dih = m_dih_anal->run(chain).score;
	m_score = m_constant + m_score_dist * m_weight_dist + m_score_dih * m_weight_dih;
}

void ScoreAa::train(const Chain &c) {
	m_dist_anal->train(c);
}

void ScoreAa::print_counts(std::ostream & stream) const {
	m_dist_anal->print_counts(stream);
}

void ScoreAa::print_freqs(std::ostream & stream) const {
	m_dist_anal->print_freqs(stream);
}

double ScoreAa::en_stacking(const Residue &r1, const Residue &r2) {
	return m_dist_anal->en_stacking(r1, r2);
}

double ScoreAa::en_pairing(const Residue &r1, const Residue &r2) {
	return m_dist_anal->en_pairing(r1, r2);

}

Score &ScoreAa::en_bp(const Residue &r1, const Residue &r2) {
	m_en_stacking = m_dist_anal->en_stacking(r1, r2);
	m_en_pairing = m_dist_anal->en_pairing(r1, r2);
	return *this;
}

double ScoreAa::en_len(const Chain &c, int beg) {
	return 0;
}

double ScoreAa::en_ang(const Chain &c, int beg) {
	return 0;
}

double ScoreAa::en_dih(const Chain &c, int beg) {
	return 0;
}

double ScoreAa::en_crash(const Residue &r1, const Residue &r2) {
	return 0;
}


}

