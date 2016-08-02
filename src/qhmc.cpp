#include "nsp.hpp"
#include <jian/nuc3d/quadruple/QHMC.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(qhmcpsb) {
    nuc3d::quadruple::QHMC<nuc3d::mc::MCpsb> qua(par);
    qua.run_job();
}

REGISTER_NSP_COMPONENT(qhmc1p) {
    nuc3d::quadruple::QHMC<nuc3d::mc::MC1p> qua(par);
    qua.run_job();
}

} // namespace jian

