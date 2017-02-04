#include "ScorePsb.hpp"

BEGIN_JN

REG_SCORER("psb", ScorePsb);

ScorePsb::ScorePsb() {
	m_cg = CG::fac_t::create("psb");
}

END_JN