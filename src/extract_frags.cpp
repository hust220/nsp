#include "nsp.hpp"
#include <jian/nuc3d/C2A.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(extract_frags) {
    C2A::extract_frags(Model(par[2]));
}

} // namespace jian

