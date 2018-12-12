#include <string>
#include "nsp.hpp"
#include "pdb.hpp"

namespace jian {

static S seq_chain(const Chain &c) {
    S seq(c.size(), 'X');
    for (int i = 0; i < seq.size(); i++) {
        seq[i] = c[i].name.back();
    }
    return seq;
}

REGISTER_NSP_COMPONENT(seq) {
    auto g = par.getv("global");
    auto && m = mol_read_to<Model>(g[1]);
    if (par.has("chain")) {
        for (auto && chain : m) {
            if (chain.name == par["chain"][0]) {
                std::cout << seq_chain(chain) << std::endl;
                return;
            }
        }
        std::cout << "This molecule has no chain named " << par["chain"][0] << std::endl;
    } else {
        std::vector<std::string> r {"A", "U", "G", "C"};
        std::vector<std::string> d {"DA", "DT", "DG", "DC"};
        S seq;
        for (auto && chain : m) {
            for (auto && res : chain) {
                seq += pdb::res_name(res.name, "abbr");
//                if (std::find(r.begin(), r.end(), res.name) != r.end()) {
//                    seq += res.name;
//                } else if (std::find(d.begin(), d.end(), res.name) != d.end()) {
//                    seq += res.name.back();
//                } else {
//                    seq += 'X';
//                }
            }
        }
        std::cout << seq << std::endl;
    }
}

}

