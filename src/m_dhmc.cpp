#include "nsp.hpp"
#include <jian/dhmc/DHMC.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(mcpsb) {
    nuc3d::DHMC<CGpsb> mcpsb;
    mcpsb.init(par);
    mcpsb.predict();
}

REGISTER_NSP_COMPONENT(mc1p) {
    nuc3d::DHMC<CG1p> mc1p;
    mc1p.init(par);
    mc1p.predict();
}

} // namespace jian
