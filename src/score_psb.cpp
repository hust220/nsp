#include "score_psb.hpp"

namespace jian {

REG_SCORER("psb", ScorePsb);

ScorePsb::ScorePsb() {
	m_cg = CG::fac_t::create("psb");
}

}
