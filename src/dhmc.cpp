#include "nsp.hpp"
#include <jian/dhmc/DHMC.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(mcpsb) {
    nuc3d::DHMC<nuc3d::mc::MCpsb> mcpsb(par);
    mcpsb.run_job();
}

REGISTER_NSP_COMPONENT(mc1p) {
    nuc3d::DHMC<nuc3d::mc::MC1p> mc1p(par);
    mc1p.run_job();
}

} // namespace jian

