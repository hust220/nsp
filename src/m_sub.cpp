#include "nsp.hpp"
#include <nsp/pdb.hpp>
#include <jian/geom.hpp>
#include <nsp/cg.hpp>

BEGIN_JN

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

    REGISTER_NSP_COMPONENT(sub) {
        std::deque<int> ls;
        Par::pars_t nums, chains;
        S in, out;
        Molecule mol;
        std::ostream stream(std::cout.rdbuf());
        std::ofstream ofile;
        tokenize_v v, w;
        int i;

        par.set(in, "p", "pdb", "s");
        mol_read(mol, in);
        if (mol.empty()) {
            throw "Empty file!";
        }

        par.set(out, "o", "out");
        if (out != "") {
            FOPEN(ofile, out);
            stream.rdbuf(ofile.rdbuf());
        }

        if (!par.has("c", "chains") && par.has("n", "num")) {
            nums = par.getv("n", "num");
            set_nums(ls, nums);
            mol = sub(mol, ls);
            //stream << sub(mol, ls) << std::endl;
        }
        else if (par.has("c", "chains")) {
            par.setv(chains, "c", "chains");
            Molecule mol_new;
            mol_new.resize(mol.size());
            for (auto && chain : chains) {
                tokenize(chain, v, ":");
                i = 0;
                for (auto && m : mol) {
                    auto it = std::find_if(m.begin(), m.end(), [&v](auto && c) {
                            return c.name == v[0];
                            });
                    if (it == m.end()) {
                        throw "The molecule has no chain named '" + v[0] + "' !";
                    }
                    else {
                        if (v.size() == 2) {
                            tokenize(v[1], w, "+");
                            set_nums(ls, w);
                            mol_new[i].push_back(sub(*it, ls));
                        }
                        else {
                            mol_new[i].push_back(*it);
                        }
                    }
                    i++;
                }
            }
            mol = mol_new;
            //stream << mol_new << std::endl;
        }
        if (par.has("cg")) {
            S cg = par.get("cg");
            CG *m_cg = CG::fac_t::create(cg);
            mol = m_cg->to_cg(mol);
            delete m_cg;
        }
        stream << mol << std::endl;

    }
}

END_JN
















