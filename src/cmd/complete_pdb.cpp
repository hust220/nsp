#include "nsp.hpp"
#include "pdb.hpp"
#include "geom.hpp"
#include "env.hpp"
#include "pdb_reader.hpp"

namespace jian {

static Chain read_template() {
    return read_model_to_chain(to_str(Env::lib(), "/RNA/pars/cg/CG2AA/templates.pdb"));
}

static geom::Superposition<Num> suppos_res(const Residue &r1, const Residue &r2) {
    List<Array<Num, 3>> ls1, ls2;
    for (auto && atom : r1) {
        auto it = std::find_if(r2.begin(), r2.end(), [&atom](auto && a){return atom.name == a.name;});
        if (it != r2.end()) {
            ls1.push_back({atom[0], atom[1], atom[2]});
            ls2.push_back({it->at(0), it->at(1), it->at(2)});
        }
    }
    Int l = size(ls1);
    Mat m1(l, 3), m2(l, 3);
    Int i = 0;
    for (auto it1 = ls1.begin(), it2 = ls2.begin(); it1 != ls1.end(); it1++, it2++) {
        for (Int j = 0; j < 3; j++) {
            m1(i, j) = it1->at(j);
            m2(i, j) = it2->at(j);
        }
        i++;
    }
    return geom::suppos(m1, m2);
}

static Bool res_is_complete(const Residue &r) {
    const auto &names = pdb::Names::instance("RNA").atoms_res.at(r.name);
    return std::all_of(names.begin(), names.end(), [&r](auto && name){return std::find_if(r.begin(), r.end(), [&name](auto && a){return name == a.name;}) != r.end();});
}

static Num rmsd_res(const Residue &r1, const Residue &r2) {
    return suppos_res(r1, r2).rmsd;
}

static void substitute_atoms(Residue &r1, const Residue r2) {
    auto sp = suppos_res(r1, r2);
    for (auto && atom : r1) {
        sp.apply(atom);
        auto it = std::find_if(r2.begin(), r2.end(), [&atom](auto && a){return atom.name == a.name;});
        if (it != r2.end()) {
            atom = *it;
        }
    }
}

static Residue complete_res(const Residue &res) {
    static Chain chain = read_template();
    Int l = size(chain);
    Vector<Num> rmsds(l);
    for (Int i = 0; i < l; i++) {
        if (chain[i].name == res.name && res_is_complete(chain[i])) {
            rmsds[i] = rmsd_res(chain[i], res);
        }
        else {
            rmsds[i] = 999;
        }
    }
    auto it = std::min_element(rmsds.begin(), rmsds.end());
    Int d = std::distance(rmsds.begin(), it);
    substitute_atoms(chain[d], res);
    return chain[d];
}

static void complete_pdb(Molecule &mol) {
    for (auto && model : mol) {
        for (auto && chain : model) {
            for (auto && res : chain) {
                if (res.name == "A" || res.name == "U" || res.name == "G" || res.name == "C") {
                    if (!res_is_complete(res)) {
                        res = complete_res(res);
                    }
                }
            }
        }
    }
}

REGISTER_NSP_COMPONENT(complete_pdb) {
    auto g = par.getv("global");

    Molecule mol;
    PdbReader reader(mol);
    reader.read(g[1]);

    complete_pdb(mol);
    JN_OUT << mol << std::endl;
}

}

