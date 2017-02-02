#include "nsp.hpp"
#include <jnbio/dhmc/DHMC.hpp>
#include <jnbio/thmc/THMC.hpp>
#include <jnbio/qhmc/QHMC.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(opt) {
	auto dhmc = DHMC::make(par);
	dhmc->predict();
}

ALIAS_NSP_COMPONENT(opt, dhmc);

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

END_JN

