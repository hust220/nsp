#include "nsp.hpp"
#include <jian/nuc3d/MCpsb.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(mcpsb) {
    nuc3d::MCpsb mcpsb(par);
    mcpsb.run();
}

} // namespace jian

