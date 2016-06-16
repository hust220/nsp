#include <iostream>
#include "nsp.hpp"
#include <jian/pdb.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(rmsd) {
    RMSD rmsd;
    std::cout << rmsd(Model(par["pdb"][0]), Model(par["pdb"][1])) << std::endl;
}

} // namespace jian

