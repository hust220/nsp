#include "nsp.hpp"
#include <jian/pdb.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(len) {
    Model m(par["pdb"][0]);
    if (par.has("chain")) {
        for (auto && chain : m) {
            if (chain.name == par["chain"][0]) {
                std::cout << chain.size() << std::endl;
                return;
            }
        }
        std::cout << "This molecule has no chain named " << par["chain"][0] << std::endl;
    } else {
        std::cout << num_residues(m) << std::endl;
    }
}

} // namespace jian

