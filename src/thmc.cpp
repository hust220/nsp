#include "nsp.hpp"
#include <jian/thmc/THMC.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(thmcpsb) {
    nuc3d::triple::THMC<nuc3d::mc::MCpsb> tri;
    tri.init(par);
    tri.run_job();
}

REGISTER_NSP_COMPONENT(thmc1p) {
    nuc3d::triple::THMC<nuc3d::mc::MC1p> tri;
    tri.init(par);
    tri.run_job();
}

} // namespace jian

