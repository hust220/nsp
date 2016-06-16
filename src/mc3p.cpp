#include "nsp.hpp"
#include <jian/nuc3d/MC3p.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(mc3p) {
    nuc3d::MC3p mc3p(par);
    mc3p.run();
}

} // namespace jian

