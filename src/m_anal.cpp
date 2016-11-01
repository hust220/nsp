#include "nsp.hpp"
#include <jian/pdb/Model.hpp>
#include <jian/geom.hpp>
#include <jian/pp.hpp>
#include <jian/utils/Debug.hpp>
#include <jian/nuc3d/ParseHelix.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(anal) {
	std::ofstream ofile;
	std::ostream stream(std::cout.rdbuf());
	Chain chain;
	int len;
    std::deque<int> ls; 
	std::string atom = "C4*";

	if (par.has("o", "out")) {
		ofile.open(par.get("o", "out"));
		stream.rdbuf(ofile.rdbuf());
	}

	chain_read_model(chain, par.get("s"));
    len = chain.size();
    if (par.has("num", "n")) {
		for (auto && s : par.getv("num", "n")) {
			ls.push_back(JN_INT(s));
		}
    }

    par.set(atom, "atom");

    if (par[2] == "dist") {
        stream << geom::distance(chain[ls[0]-1][atom], chain[ls[1]-1][atom]) << std::endl;
    } else if (par[2] == "dist_atom") {
        stream << geom::distance(chain[JN_INT(par["atom1"][0])-1][par["atom1"][1]],
			                     chain[JN_INT(par["atom2"][0])-1][par["atom2"][1]]) << std::endl;
    } else if (par[2] == "ang") {
        if (par.has("all")) {
            std::deque<Residue> dq;
            for (int i = 0; i < len; i++) {
                if ((!dq.empty()) && geom::distance(dq.back()["O3*"], chain[i]["C5*"]) > 4) {
                    dq.clear();
                }
                dq.push_back(chain[i]);
                if (dq.size() == 3) {
                    stream << geom::angle(dq[0][atom], dq[1][atom], dq[2][atom]) << std::endl;
                    dq.pop_front();
                }
            }
        } else {
            stream << geom::angle(chain[ls[0]-1][atom], chain[ls[1]-1][atom], chain[ls[2]-1][atom]) << std::endl;
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
                    stream << geom::dihedral(dq[0][atom], dq[1][atom], dq[2][atom], dq[3][atom]) << std::endl;
                    dq.pop_front();
                }
            }
        } else {
            stream << geom::dihedral(chain[ls[0]-1][atom], chain[ls[1]-1][atom],
				                     chain[ls[2]-1][atom], chain[ls[3]-1][atom]) << std::endl;
        }
    } else if (par[2] == "chir") {
        stream << geom::chirality(chain[ls[0]-1][atom], chain[ls[1]-1][atom], chain[ls[2]-1][atom], chain[ls[3]-1][atom]) << std::endl;
    }
}

REGISTER_NSP_COMPONENT(dihs) {
    std::string pdb = par.get("s");
    std::deque<Atom> atoms;
    Chain chain = read_model_to_chain(pdb);
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
















