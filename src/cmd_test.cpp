#include "nsp.hpp"
#include "pdb.hpp"
#include "geom.hpp"

namespace jian {

namespace {

template<typename NUMS>
void set_nums(std::deque<int> &ls, NUMS &&nums) {
    tokenize_v v;
    int beg, end;

    for (auto && s : nums) {
        tokenize(s, v, "-");
        if (v.size() == 1) {
            ls.push_back(JN_INT(s) - 1);
        }
        else if (v.size() == 2) {
            beg = JN_INT(v[0]) - 1;
            end = JN_INT(v[1]) - 1;
            for (int i = beg; i <= end; i++) {
                ls.push_back(i);
            }
        }
    }
}

REGISTER_NSP_COMPONENT(test) {
    auto g = par.getv("global");
    auto && m = mol_read_to<Molecule>(g[1]);
    Str a = g[2];
    Str b = g[3];

    auto center = [](auto atoms, auto s) {
        tokenize_v v, w;

        tokenize(s, v, "+");
        Deque<Int> ls;
        if (v.size() == 1) {
            ls.push_back(JN_INT(v[0])-1);
        }
        else {
            for (auto && i : v) {
                tokenize(i, w, "-");
                for (int j = JN_INT(w[0])-1; j < JN_INT(w[1]); j++) {
                    ls.push_back(j);
                }
            }
        }

        //for (auto && n : ls) JN_OUT << n << ' '; JN_OUT << Endl;

        Vec c = Vec::Zero(3);
        for (int i = 0; i < 3; i++) {
            for (auto && n : ls) {
                c[i] += atoms[n][i];
            }
            c[i] /= double(ls.size());
        }

        return c;
    };

    for (auto && model : m) {
        auto rs = model.residues();
        Deque<Atom> atoms;
        for (const Residue & r : rs) {
            atoms.push_back(r["CA"]);
        }


        auto v1 = center(atoms, a);
        //JN_OUT << "v1: " << v1 << Endl;
        auto v2 = center(atoms, b);
        //JN_OUT << "v2: " << v2 << Endl;

        JN_OUT << geom::distance(v1, v2) << Endl;
    }
}

}

}

