#include "nsp.hpp"
#include <jian/nuc3d/triple/THMC.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(thmcpsb) {
    nuc3d::triple::THMC<nuc3d::mc::MCpsb> tri(par);
    tri.run_job();
}

REGISTER_NSP_COMPONENT(thmc1p) {
    nuc3d::triple::THMC<nuc3d::mc::MC1p> tri(par);
    tri.run_job();
}

} // namespace jian

