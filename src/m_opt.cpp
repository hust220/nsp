#include "nsp.hpp"
#include <jian/dhmc/DHMC.hpp>
#include <jian/thmc/THMC.hpp>
#include <jian/qhmc/QHMC.hpp>

namespace jian {

	REGISTER_NSP_COMPONENT(dhmc) {
		DHMC dhmc;
		dhmc.init(par);
		dhmc.predict();
	}

	REGISTER_NSP_COMPONENT(thmc) {
		nuc3d::triple::THMC tri;
		tri.init(par);
		tri.predict();
	}

	REGISTER_NSP_COMPONENT(qhmc) {
		qhmc::QHMC qua;
		qua.init(par);
		qua.predict();
	}

} // namespace jian

