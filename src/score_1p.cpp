#include "score_1p.hpp"

namespace jian {

REG_SCORER("1p", Score1p);

Score1p::Score1p() {
	m_cg = CG::fac_t::create("1p");
}

}
