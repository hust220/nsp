#include <string>
#include "../utils/Env.hpp"
#include "../utils/Par.hpp"
#include "../pdb.hpp"
#include "DistAnal.hpp"
#include "DihAnal.hpp"
#include "Scoring.hpp"

namespace jian {

Scoring::Scoring() {
    std::string lib = Env::lib() + "/RNA/pars/scoring";

    m_par_dist = lib + "/dist.frequencies";
    m_par_dih  = lib + "/par_dih";

    m_dist_anal = new DistAnal(m_bin_dist, m_cutoff);
    m_dist_anal->read_parm(m_par_dist);

    m_dih_anal = new DihAnal(m_bin_dih);
    m_dih_anal->read_parm(m_par_dih);
}

Scoring::~Scoring() {
    if (m_dist_anal != NULL) delete m_dist_anal;
    if (m_dih_anal != NULL) delete m_dih_anal;
}

Scoring &Scoring::run(const Chain &chain) {
    m_score_dist = m_dist_anal->run(chain).score;
    m_score_dih = m_dih_anal->run(chain).score;
    m_score = m_constant + m_score_dist * m_weight_dist + m_score_dih * m_weight_dih;
    return *this;
}

Scoring & train(std::string filename) {
	return *this;
}

} // namespace jian

