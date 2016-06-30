#include "nsp.hpp"
#include <jian/pdb/Model.hpp>
#include <jian/geom.hpp>
#include <jian/pp.hpp>
#include <jian/utils/Debug.hpp>
#include <jian/nuc3d/ParseHelix.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(anal_mol) {
    auto &&chain = residues_from_file(par["pdb"][0]);
    int len = chain.size();
    std::deque<int> ls; 
    if (par.has("num")) {
        EACH(s, par["num"], ls.push_back(JN_INT(s)));
    }
    std::string atom = "C4*";
    par.set(atom, "atom");
    if (par[2] == "dist") {
        Debug::println(geom::distance(chain[ls[0]-1][atom], chain[ls[1]-1][atom]));
    } else if (par[2] == "dist_atom") {
        Debug::println(geom::distance(chain[JN_INT(par["atom1"][0])-1][par["atom1"][1]], chain[JN_INT(par["atom2"][0])-1][par["atom2"][1]]));
    } else if (par[2] == "ang") {
//        Debug::println(geom::angle(chain[ls[0]-1][atom], chain[ls[1]-1][atom], chain[ls[2]-1][atom]));
        if (par.has("all")) {
            std::deque<Residue> dq;
            for (int i = 0; i < len; i++) {
                if ((!dq.empty()) && geom::distance(dq.back()["O3*"], chain[i]["C5*"]) > 4) {
                    dq.clear();
                }
                dq.push_back(chain[i]);
                if (dq.size() == 3) {
                    std::cout << geom::angle(dq[0][atom], dq[1][atom], dq[2][atom]) << std::endl;
                    dq.pop_front();
                }
            }
        } else {
            Debug::println(geom::angle(chain[ls[0]-1][atom], chain[ls[1]-1][atom], chain[ls[2]-1][atom]));
        }
    } else if (par[2] == "dih") {
        if (par.has("all")) {
            std::deque<Residue> dq;
            for (int i = 0; i < len; i++) {
                if ((!dq.empty()) && geom::distance(dq.back()["O3*"], chain[i]["C5*"]) > 4) {
                    dq.clear();
                }
                dq.push_back(chain[i]);
                if (dq.size() == 4) {
                    std::cout << geom::dihedral(dq[0][atom], dq[1][atom], dq[2][atom], dq[3][atom]) << std::endl;
                    dq.pop_front();
                }
            }
        } else {
            Debug::println(geom::dihedral(chain[ls[0]-1][atom], chain[ls[1]-1][atom], chain[ls[2]-1][atom], chain[ls[3]-1][atom]));
        }
    } else if (par[2] == "chir") {
        Debug::println(geom::chirality(chain[ls[0]-1][atom], chain[ls[1]-1][atom], chain[ls[2]-1][atom], chain[ls[3]-1][atom]));
    }
}

//REGISTER_NSP_COMPONENT(num_residues) {
//    std::cout << num_residues(Model(par[2])) << std::endl;
//}

REGISTER_NSP_COMPONENT(dihs) {
    std::string pdb = par["pdb"][0];
    std::deque<Atom> atoms;
    Chain chain = residues_from_file(pdb);
    for (auto && r : chain) {
        for (auto && a : r) {
            if (a.name == "C5*" || a.name == "O3*") {
                atoms.push_back(a);
            }
        }
    }
    for (int i = 0; i < atoms.size() - 3; i++) {
        std::cout << atoms[i].name << ' ' << geom::dihedral(atoms[i], atoms[i+1], atoms[i+2], atoms[i+3]) << std::endl;
    }
}

REGISTER_NSP_COMPONENT(dihs_std_helix) {
    dihs_std_helix();
}

} // namespace jian
















