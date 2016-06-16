#include "nsp.hpp"
#include <jian/nuc3d/triple/TriPred.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(tripred) {
    TriPred tri(par);
    tri.predict();
}

} // namespace jian

