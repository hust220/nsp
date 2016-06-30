#include "CG1p.hpp"

namespace jian {

int CG1p::size_res = 1;

Residue CG1p::res(const Residue &r) {
    Residue res;
    res.name = r.name;
    for (auto && atom : r) {
        if (atom.name == "C4*") {
            res.push_back(atom);
        }
    }
    return res;
}

Chain CG1p::chain(const Chain &chain) {
    Chain c;
    c.name = chain.name;
    c.model_name = chain.model_name;
    for (auto && r : chain) {
        c.push_back(res(r));
    }
    return c;
}

} // namespace jian
