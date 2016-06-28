#include "nsp.hpp"
#include <jian/nuc3d/quadruple/QuaPred.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(qua_mcpsb) {
    nuc3d::quadruple::QuaPred<nuc3d::mc::MCpsb> qua(par);
    qua.run();
}

REGISTER_NSP_COMPONENT(qua_mc1p) {
    nuc3d::quadruple::QuaPred<nuc3d::mc::MC1p> qua(par);
    qua.run();
}

} // namespace jian

