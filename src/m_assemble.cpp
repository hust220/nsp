#include "nsp.hpp"
#include <jian/nuc3d/Assemble.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(assemble) {
    nuc3d::Assemble ass(par);
    ass.predict();
}

} // namespace jian

