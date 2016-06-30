#include "nsp.hpp"
#include <jian/pdb.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(chain) {
    Model m(par["pdb"][0]);
    for (auto && chain : m) {
        if (chain.name == par["chain"][0]) {
            std::cout << chain << std::endl;
            return;
        }
    }
    std::cout << "This molecule has no chain named " << par["chain"][0] << std::endl;
}

} // namespace jian
