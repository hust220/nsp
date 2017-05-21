#include "tsdna_set_mvels.hpp"

BEGIN_JN

namespace tsdna {

    void set_related_residues(THMC &m) {
        m.d_mc_related_residues.resize(size(m._seq));
        for (auto && module : m.d_modules) {
            if (module->type() != "helix") {
                for (auto && frag : module->d_frags) {
                    for (auto && i : frag) {
                        m.d_mc_related_residues[i] = std::make_shared<std::deque<int>>();
                        m.d_mc_related_residues[i]->push_back(i);
                    }
                }
            }
        }
        for (auto && module : m.d_modules) {
            if (module->type() == "helix") {
                auto p = std::make_shared<Deque<Int>>();
                for (auto && frag : module->d_frags) {
                    for (auto && i : frag) {
                        m.d_mc_related_residues[i] = p;
                        p->push_back(i);
                    }
                }
            }
        }
    }

    void set_unrelated_residues(THMC &m) {
        int len = size(m._seq);
        auto &r = m.d_mc_related_residues;
        m.d_mc_unrelated_residues.resize(len);
        for (int i = 0; i < len; i++) {
            m.d_mc_unrelated_residues[i] = std::make_shared<Deque<Int>>();
            for (int j = 0; j < len; j++) {
                if (std::none_of(r[i]->begin(), r[i]->end(), [&j](auto && n){return n == j;})) {
                    m.d_mc_unrelated_residues[i]->push_back(j);
                }
            }
        }
    }

    void set_mvels_helices(THMC &m) {
        for (auto && module : m.d_modules) {
            if (module->type() == "helix") {
                MvEl *el = new MvEl(MVEL_HL);
                for (auto && frag : module->d_frags) {
                    el->add_frag(frag.front(), frag.back());
                }
                m.m_mvels.push_back(el);
            }
        }
    }

    Bool is_fixed(const THMC &m, Int i) {
        return m._ss[i] != '.' && m._ss[i] != '0';
    }

    void set_mvels_frags(THMC &m) {
        Int i, j, l = size(m._seq);
        for (i = 0; i + 2 < l; i++) {
            if (!is_fixed(m, i) || !is_fixed(m, i+1) || !is_fixed(m, i+2)) {
                m.m_mvels.push_back(new MvEl(i, i + 2, MVEL_FG));
            }
        }
    }

    void thmc_set_mvels(THMC &m) {
        set_mvels_helices(m);
        set_mvels_frags(m);
    }

}

END_JN

