#include "nsp.hpp"
#include <jian/pdb.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(rna) {
    Model m(par["pdb"][0]);
    std::cout << m << std::endl;
}

} // namespace jian

