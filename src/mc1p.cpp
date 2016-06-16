#include "nsp.hpp"
#include <jian/nuc3d/MC1p.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(mc1p) {
    nuc3d::MC1p mc1p(par);
    mc1p.run();
}

} // namespace jian

