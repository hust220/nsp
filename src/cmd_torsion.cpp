#include "nsp.hpp"
#include "pdb.hpp"
#include "geom.hpp"
#include "rss_get_ss.hpp"

BEGIN_JN

namespace {

    Bool is_close(const Residue &r1, const Residue &r2) {
        return geom::distance(r1["O3*"], r2["C5*"]) < 4;
    }

    void add_res(Deque<Residue> &ls, Deque<Char> &dq_ss, Char ss, const Residue &r, int n) {
        Num chi;
        ls.push_back(r);
        dq_ss.push_back(ss);
        if (size(ls) == n) {
            Int b = n/2, a = b - 1, c = b + 1;
            Num eta = geom::dihedral(ls[a]["C4*"], ls[b]["P"], ls[b]["C4*"], ls[c]["P"]);
            Num theta = geom::dihedral(ls[b]["P"], ls[b]["C4*"], ls[c]["P"], ls[c]["C4*"]);
            if (std::find_if(ls[b].begin(), ls[b].end(), [](auto &&atom){return atom.name == "N9";})!=ls[b].end())
                chi = geom::dihedral(ls[b]["C2*"], ls[b]["C1*"], ls[b]["N9"], ls[b]["C4"]);
            else
                chi = geom::dihedral(ls[b]["C2*"], ls[b]["C1*"], ls[b]["N1"], ls[b]["C2"]);
            for (auto && i : ls) { JN_OUT << i.name; }
            JN_OUT << ' ';
            for (auto && i : dq_ss) { JN_OUT << i; }
            JN_OUT << ' ' << eta << ' ' << theta << ' ' << chi << std::endl;
            ls.pop_front();
            dq_ss.pop_front();
        }
    }

    REGISTER_NSP_COMPONENT(torsion) {
        auto global = par.getv("global");
        Str filename = global[1];
        Int n = JN_INT(global[2]);
        if (n % 2 == 0 || n < 3) throw "Illegal length of fragment!";

        Chain c;
        chain_read_model(c, filename);
        Str ss = get_ss(c);

        Deque<Residue> ls;
        Deque<Char> dq_ss;
        Int i = 0;
        for (auto && r : c) {
            if (ls.empty() || is_close(ls.back(), r)) {
                add_res(ls, dq_ss, ss[i], r, n);
            }
            else {
                ls.clear();
                dq_ss.clear();
            }
            i++;
        }
    }

}

END_JN

