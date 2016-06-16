#include "nsp.hpp"
#include <jian/pdb.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(rna) {
    Model m(par[2]);
    std::cout << m << std::endl;
}

} // namespace jian

