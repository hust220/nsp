#include "nsp.hpp"
#include <jian/nuc3d/BuildLoop.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(build_loop) {
    int n = JN_INT(par["num"][0]);
    BuildLoop build_loop;
    for (int i = 0; i < n; i++) {
        mol_write(build_loop(par["seq"][0], par["ss"][0]), par["name"][0]+'-'+JN_STR(i+1)+".pdb");
    }
}

} // namespace jian

