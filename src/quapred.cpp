#include "nsp.hpp"
#include <jian/nuc3d/quadruple/QuaPred.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(quapred) {
    QuaPred tri(par);
    tri.predict();
}

} // namespace jian

