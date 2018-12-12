#include <vector>
#include <string>
#include <set>
#include <array>
#include <algorithm>
#include "cg_aa.hpp"

namespace jian {
	REG_CG("aa", CGaa);

	CGaa::CGaa() {
		m_cg = "aa";
	}

	int CGaa::res_size() const {
		return -1;
	}

	Residue CGaa::to_cg(const Residue &r) const {
		return r;
	}

}
