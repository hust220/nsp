#include "nsp.hpp"
#include "pdb.hpp"

BEGIN_JN

REGISTER_NSP_COMPONENT(len) {
    auto g = par.getv("global");
    auto && m = mol_read_to<Model>(g[1]);
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

END_JN

