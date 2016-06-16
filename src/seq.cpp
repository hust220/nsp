#include <string>
#include "nsp.hpp"
#include <jian/pdb.hpp>

namespace jian {

static std::string seq_chain(const Chain &c) {
    std::string seq(c.size(), 'X');
    for (int i = 0; i < seq.size(); i++) {
        seq[i] = c[i].name.back();
    }
    return seq;
}

REGISTER_NSP_COMPONENT(seq) {
    Model m(par["pdb"][0]);
    if (par.has("chain")) {
        for (auto && chain : m) {
            if (chain.name == par["chain"][0]) {
                std::cout << seq_chain(chain) << std::endl;
                return;
            }
        }
        std::cout << "This molecule has no chain named " << par["chain"][0] << std::endl;
    } else {
        std::cout << seq(m) << std::endl;
    }
}

} // namespace jian

