#include "nsp.hpp"
#include <jian/geom.hpp>
#include <jian/pdb.hpp>
#include <jian/scoring/Scoring.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(score) {
    std::string s_file = par["s"][0];

    Chain chain;
    chain_read_model(chain, s_file);

    Scoring scoring;
    scoring.run(chain);
    std::cout << scoring.m_score_dih << "(dih) "
              << scoring.m_score_dist << "(dist) "
              << scoring.m_score << "(total)"
              << std::endl;
}

} // namespace jian
















