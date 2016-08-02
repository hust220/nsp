#include "nsp.hpp"
#include <jian/nuc3d/split.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(split) {
    nuc3d::split(par);
}

} // namespace jian

