#include "nsp.hpp"
#include <nsp/pdb/Model.hpp>
#include <jian/geom.hpp>
#include <jian/pp.hpp>
#include <jian/utils/Debug.hpp>
#include <nsp/nuc3d/ParseHelix.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(anal) {
	Chain chain;
	int len;
    std::deque<int> ls; 
	S atom = "C4*";

	chain_read_model(chain, par.get("s"));
    len = chain.size();
    if (par.has("num", "n")) {
		for (auto && s : par.getv("num", "n")) {
			ls.push_back(JN_INT(s));
		}
    }

    par.set(atom, "atom");

	std::deque<int> num_res;
	std::deque<std::string> name_atom;
	tokenize_v v;
	if (par.has("atoms")) {
		for (auto && s : par.getv("atoms")) {
			tokenize(s, v, ":");
			num_res.push_back(JN_INT(v[0]) - 1);
			name_atom.push_back(v[1]);
		}
	}

    if (par[2] == "dist") {
        JN_OUT << geom::distance(chain[ls[0]-1][atom], chain[ls[1]-1][atom]) << std::endl;
	}
	else if (par[2] == "dist_atom") {
		JN_OUT <<
			chain[num_res[0]].name << num_res[0] + 1 << ' ' <<
			chain[num_res[1]].name << num_res[1] + 1 << ' ' <<
			geom::distance(
				chain[num_res[0]][name_atom[0]],
				chain[num_res[1]][name_atom[1]]) <<
			std::endl;
	}
	else if (par[2] == "ang_atom") {
		JN_OUT <<
			chain[num_res[0]].name << num_res[0] + 1 << ' ' <<
			chain[num_res[1]].name << num_res[1] + 1 << ' ' <<
			chain[num_res[2]].name << num_res[2] + 1 << ' ' <<
			geom::angle(
				chain[num_res[0]][name_atom[0]],
				chain[num_res[1]][name_atom[1]],
				chain[num_res[2]][name_atom[2]]) <<
			std::endl;
	}
	else if (par[2] == "dih_atom") {
		JN_OUT <<
			chain[num_res[0]].name << num_res[0] + 1 << ' ' <<
			chain[num_res[1]].name << num_res[1] + 1 << ' ' <<
			chain[num_res[2]].name << num_res[2] + 1 << ' ' <<
			chain[num_res[3]].name << num_res[3] + 1 << ' ' <<
			geom::dihedral(
				chain[num_res[0]][name_atom[0]],
				chain[num_res[1]][name_atom[1]],				
				chain[num_res[2]][name_atom[2]],
				chain[num_res[3]][name_atom[3]]) <<
			std::endl;
	}
	else if (par[2] == "ang") {
        if (par.has("all")) {
            std::deque<Residue> dq;
            for (int i = 0; i < len; i++) {
                if ((!dq.empty()) && geom::distance(dq.back()["O3*"], chain[i]["C5*"]) > 4) {
                    dq.clear();
                }
                dq.push_back(chain[i]);
                if (dq.size() == 3) {
					JN_OUT << geom::angle(dq[0][atom], dq[1][atom], dq[2][atom]) << std::endl;
                    dq.pop_front();
                }
            }
        } else {
			JN_OUT << geom::angle(chain[ls[0]-1][atom], chain[ls[1]-1][atom], chain[ls[2]-1][atom]) << std::endl;
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
					JN_OUT << geom::dihedral(dq[0][atom], dq[1][atom], dq[2][atom], dq[3][atom]) << std::endl;
                    dq.pop_front();
                }
            }
        } else {
			JN_OUT << geom::dihedral(chain[ls[0]-1][atom], chain[ls[1]-1][atom],
				                     chain[ls[2]-1][atom], chain[ls[3]-1][atom]) << std::endl;
        }
    } else if (par[2] == "chir") {
		JN_OUT << geom::chirality(chain[ls[0]-1][atom], chain[ls[1]-1][atom], chain[ls[2]-1][atom], chain[ls[3]-1][atom]) << std::endl;
    }
}

REGISTER_NSP_COMPONENT(dihs) {
    S pdb = par.get("s");
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

END_JN
















