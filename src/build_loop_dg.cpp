#include "nsp.hpp"
#include <jian/nuc3d/BuildLoopDG.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(build_loop_dg) {
    BuildLoopDG build_loop_dg;
    build_loop_dg.init(par["seq"][0], par["ss"][0]);
//    write_pdb(build_loop_dg(), "3drna-dg.pdb");
    std::cout << build_loop_dg() << std::endl;
}

} // namespace jian

