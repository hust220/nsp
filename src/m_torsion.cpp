#include "nsp.hpp"
#include <nsp/pdb.hpp>
#include <jian/geom.hpp>

BEGIN_JN

namespace {

    Bool is_close(const Residue &r1, const Residue &r2) {
        return geom::distance(r1["O3*"], r2["C5*"]) < 4;
    }

    void add_res(Deque<Residue> &ls, const Residue &r, int n) {
        ls.push_back(r);
        if (size(ls) == n) {
            Int b = n/2, a = b - 1, c = b + 1;
            Num eta = geom::dihedral(ls[a]["C4*"], ls[b]["P"], ls[b]["C4*"], ls[c]["P"]);
            Num theta = geom::dihedral(ls[b]["P"], ls[b]["C4*"], ls[c]["P"], ls[c]["C4*"]);
            for (auto && i : ls) {
                JN_OUT << i.name;
            }
            JN_OUT << ' ' << eta << ' ' << theta << std::endl;
            ls.pop_front();
        }
    }

    REGISTER_NSP_COMPONENT(torsion) {
        auto global = par.getv("global");
        Str filename = global[1];
        Int n = JN_INT(global[2]);
        if (n % 2 == 0 || n < 3) throw "Illegal length of fragment!";

        Chain c;
        mol_read(c, filename);

        Deque<Residue> ls;
        for (auto && r : c) {
            if (ls.empty() || is_close(ls.back(), r)) {
                add_res(ls, r, n);
            }
        }
    }

}

END_JN

