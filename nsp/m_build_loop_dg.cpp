#include "nsp.hpp"
#include <jnbio/nuc3d/BuildLoopDG.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(build_loop_dg) {
    BuildLoopDG build_loop_dg;
    build_loop_dg.init(par["seq"][0], par["ss"][0]);
//    write_pdb(build_loop_dg(), "3drna-dg.pdb");
    std::cout << build_loop_dg() << std::endl;
}

END_JN

