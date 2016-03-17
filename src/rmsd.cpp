#include "nsp.hpp"
#include <jian/pdb/RMSD.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(rmsd) {
    RMSD rmsd;
    std::cout << rmsd(Model(par[2]), Model(par[3])) << std::endl;
}

} // namespace jian

