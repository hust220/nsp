#include "nsp.hpp"
#include <jian/pdb/Model.hpp>
#include <jian/geom.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(seq) {
    std::cout << seq(Model(par[2])) << std::endl;
}

} // namespace jian

