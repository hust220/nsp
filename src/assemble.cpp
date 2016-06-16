#include "nsp.hpp"
#include <jian/nuc3d/Ass.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(assemble) {
    nuc3d::Assemble ass(par);
    ass.run();
}

} // namespace jian

