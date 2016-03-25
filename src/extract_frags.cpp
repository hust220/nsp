#include "nsp.hpp"
#include <jian/nuc3d/C2A.hpp>
#include <jian/nuc3d/CG2AA.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(cg2aa) {
    std::vector<std::string> _suppos_atoms {"C5*", "O3*", "C1*"};
    auto &&chain = residues_from_file(par[2]);
    int len = chain.size();
    int num_atoms = len * _suppos_atoms.size();
    MatrixXd c(num_atoms, 3);
    int num_atom = 0;
    for (int i = 0; i < chain.size(); i++) {
        chain[i].sort();
        for (auto && atom : chain[i]) {
            if (std::find(_suppos_atoms.begin(), _suppos_atoms.end(), atom.name) != _suppos_atoms.end()) {
                for (int j = 0; j < 3; j++) {
                    c(num_atom, j) = atom[j];
                }
                num_atom++;
            }
        }
    }
    auto && new_chain = CG2AA::cg2aa(c, std::array<int, 2>{0, num_atoms-1});
    residues_to_file(new_chain, par[3]);
}

REGISTER_NSP_COMPONENT(c2a_extract_frags) {
    C2A::extract_frags(Model(par[2]));
}

REGISTER_NSP_COMPONENT(cg2aa_extract_frags) {
    CG2AA::extract_frags(Model(par[2]));
}

} // namespace jian

