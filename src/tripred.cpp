#include "nsp.hpp"
#include <jian/nuc3d/triple/TriPred.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(tri_mcpsb) {
    nuc3d::triple::TriPred<nuc3d::mc::MCpsb> tri(par);
    tri.run();
}

REGISTER_NSP_COMPONENT(tri_mc1p) {
    nuc3d::triple::TriPred<nuc3d::mc::MC1p> tri(par);
    tri.run();
}

} // namespace jian

