#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/geom.hpp>
#include <jian/nuc2d/SSTree.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(rmsd) {
    RMSD rmsd;
    std::cout << rmsd(Model(par[2]), Model(par[3])) << std::endl;
}

REGISTER_NSP_COMPONENT(seq) {
    std::cout << seq(Model(par[2])) << std::endl;
}

REGISTER_NSP_COMPONENT(ss_tree) {
    try {
        SSTree ss_tree;
        ss_tree.make(par["seq"][0], par["ss"][0]);
        ss_tree.head->print_tree();
    } catch(const Error &e) {
        std::cout << e.what() << std::endl;
    }
}

REGISTER_NSP_COMPONENT(view_mol_constraints) {
    auto &&chain = residues_from_file(par["m"][0]);
    int i, j;
    auto min_distance = [](auto &&r1, auto &&r2){
        double min = 999999;
        double d;
        for (auto && atom1 : r1) {
            for (auto && atom2 : r2) {
                d = geom::distance(atom1, atom2);
                if (d < min) {
                    min = d;
                }
            }
        }
        return min;
    };
    EACH_SPLIT_LINE(par["f"][0].c_str(), " ",
        if (F.size() >= 2) {
            i = JN_INT(F[0]) - 1;
            j = JN_INT(F[1]) - 1;
            std::cout << i << ' ' << j << ' ' << min_distance(chain[i], chain[j]) << std::endl;
        }
    );
}

} // namespace jian

