#include "Score1p.hpp"

BEGIN_JN

REG_SCORER("1p", Score1p);

Score1p::Score1p() {
	m_cg = CG::fac_t::create("1p");
}

END_JN