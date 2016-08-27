#pragma once

#include <iostream>
#include <set>
#include <memory>
#include <sstream>
#include "../mc.hpp"
#include "../utils/file.hpp"
#include "../scoring/score_psb.hpp"
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../nuc2d.hpp"
#include "../cg.hpp"
#include "../pp.hpp"
#include "BuildHelix.hpp"
#include "transform.hpp"
#include "TemplRec.hpp"
#include "JobPredict3D.hpp"
#include "mc/MC1p.hpp"
#include "mc/MCpsb.hpp"

namespace jian {
namespace nuc3d {

using fixed_ranges_t = std::list<std::array<int, 2>>;

template<typename MC_T>
class DHMC : public MC_T {
public:    
    using mc_t = MC_T;
    using range_t = std::array<int, 4>;

    std::deque<std::shared_ptr<SSTree>> m_trees;
    std::deque<range_t *> m_range;
    range_t _moved_residues;
    fixed_ranges_t _fixed_ranges;

    template<typename T>
    std::string partial_ss(std::string ss, T &&pair) {
        for (auto && c : ss) {
            if (c == pair.first) {
                c = '(';
            } else if (c == pair.second) {
                c = ')';
            } else if (c != '&') {
                c = '.';
            }
        }
        return ss;
    }

    template<typename T>
    auto make_hp_range(T &&l) {
        return make_range(l->s.head->res1.num - 1, l->s.head->res2.num - 1, l->s.head->res1.num - 1, l->s.head->res2.num - 1 );
    }

    template<typename T>
    auto make_il_range(T &&l) {
        return make_range(l->s.head->res1.num - 1, l->son->s.head->res1.num - 2, l->son->s.head->res2.num, l->s.head->res2.num - 1 );
    }

    template<typename T>
    auto make_helix_range(T &&h) {
        auto p = make_range(h.head->res1.num - 1, 0, 0, h.head->res2.num - 1);
        HELIX_EACH(h,
            if (BP->next == NULL) {
                (*p)[1] = BP->res1.num - 1;
                (*p)[2] = BP->res2.num - 1;
            }
        );
        return p;
    }

    DHMC(const Par &par) : mc_t(par) {
        LOG << "# Set 2D trees" << std::endl;
        set_trees();

        LOG << "# Set ranges" << std::endl;
        set_ranges();

        LOG << "# Print ranges" << std::endl;
        print_ranges();
    }

    ~DHMC() {
        for (auto && r : m_range) {
            delete r;
        }
    }

    void print_ranges() {
        for (auto && range : m_range) {
            for (auto && i : *range) {
                LOG << i << ' ';
            }
            LOG << std::endl;
        }
    }

    void set_trees() {
        m_trees.push_back(std::make_shared<SSTree>());
        m_trees.back()->make_b(mc_t::_seq, mc_t::_ss, 2);
        auto & keys = NucSS::instance().paired_keys;
        for (auto it = keys.begin()+1; it != keys.end(); it++) {
            auto ss = partial_ss(mc_t::_ss, *it);
            if (std::any_of(ss.begin(), ss.end(), [](auto &&c){return c != '.' && c != '&';})) {
                m_trees.push_back(std::make_shared<SSTree>());
                m_trees.back()->make_b(mc_t::_seq, ss, 1);
            } else {
                break;
            }
        }
    }

    void translate_pseudo_knots_helix(auto & m, auto && nums) {
        int n = mc_t::cg_t::size_res * nums.size() / 2;
        Mat x(n, 3), y(n, 3);
        int i = 0;
        int l = 0;
        for (auto && j : nums) {
            if (i < nums.size() / 2) {
                auto && res1 = mc_t::cg_t::res(m[0][i]);
                auto && res2 = mc_t::cg_t::res(mc_t::_pred_chain[j]);
                for (int k = 0; k < mc_t::cg_t::size_res; k++) {
                    for (int t = 0; t < 3; t++) {
                        x(l, t) = res1[k][t];
                        y(l, t) = res2[k][t];
                    }
                    l++;
                }
            }
            i++;
        }
        auto sp = geom::suppos(x, y);
        INIT_SUPPOS(sp);
        for (auto && res : m[0]) {
            for (auto && atom : res) {
                APPLY_SUPPOS(atom, sp);
            }
        }
    }

    void set_pseudo_knots_helix (auto &&h){
        auto && seq = h.seq();
        auto && m = build_helix(seq);
        auto && nums = h.nums();

        assert(nums.size() >= 2 && nums.size() % 2 == 0);
        translate_pseudo_knots_helix(m, nums);

        int i = 0;
        for (auto && n : nums) {
            mc_t::_pred_chain[n] = std::move(m[0][i]);
            i++;
        }
    }

    void set_pseudo_knots() {
        auto it = m_trees.begin();

        LOOP_TRAVERSE((*it)->head(), 
            if (L->has_helix() && L->s.len() == 1) {
                set_pseudo_knots_helix(L->s);
            }
        );

        for (it = m_trees.begin() + 1; it != m_trees.end(); it++) {
            LOOP_TRAVERSE((*it)->head(), 
                if (L->has_helix()) {
                    set_pseudo_knots_helix(L->s);
                }
            );
        }
    }

    auto make_range(int a, int b, int c,  int d) {
//        std::shared_ptr<range_t> p = std::make_shared<range_t>();
        range_t *p = new range_t;
        (*p) = {a, b, c, d};
        return p;
    }

    void set_ranges() {
        std::vector<int> v(mc_t::_seq.size(), 0);

        auto is_hp = [&](loop *l) {
            if (l->has_loop() && l->has_helix() && l->num_sons() == 0) {
                int flag = 0, n = 0;
                LOOP_EACH(l,
                    char &c = mc_t::_ss[RES->num - 1];
                    if (c != '.'  && c != '(' && c != ')') {
                        flag = 1;
                        break;
                    } else {
                        n++;
                    }
                );
                if (flag == 0 && n <= 14) {
                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        };

        auto is_il = [&](loop *l) {
            if (l->has_loop() && l->has_helix() && l->num_sons() == 1) {
                int flag = 0, n = 0;
                LOOP_EACH(l,
                    char &c = mc_t::_ss[RES->num - 1];
                    if (c != '.'  && c != '(' && c != ')') {
                        flag = 1;
                        break;
                    } else {
                        n++;
                    }
                );
                if (flag == 0 && n <= 20) {
                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        };

        auto update_range = [&](auto &&range) {
            for (int i = range->at(0); i <= range->at(1); i++) { v[i] = 1; }
            for (int i = range->at(2); i <= range->at(3); i++) { v[i] = 1; }
            m_range.push_back(range);
        };

        auto set_res_module_types_ss = [&](loop *l, bool is_first){
            LOOP_TRAVERSE(l,
                if (!mc_t::m_sample_hairpin && is_first && is_hp(L)) {
                    update_range(make_hp_range(L));
                } else if (!mc_t::m_sample_il && is_first && is_il(L)) {
                    update_range(make_il_range(L));
                } else if (L->has_helix()) {
                    update_range(make_helix_range(L->s));
                }
            );
        };

        for (auto it = m_trees.begin(); it != m_trees.end(); it++) {
            set_res_module_types_ss((*it)->head(), it == m_trees.begin());
        }

        for (int i = 0; i < mc_t::_seq.size(); i++) {
            if (v[i] == 0) {
                m_range.push_back(make_range(i, i, i, i));
            }
        }

        merge_ranges();
        set_fixed_ranges();
    }

    bool in_fixed_ranges(const range_t &r) {
        for (auto && range : _fixed_ranges) {
            if ((range[0] >= r[1] && range[1] >= r[0]) || 
                (range[0] >= r[3] && range[1] >= r[2])) {
                return true;
            }
        }
        return false;
    }

    void set_fixed_ranges() {
        while (true) {
            bool flag = false;
            for (auto && r : m_range) {
                auto & range = *r;
                if (in_fixed_ranges(range)) {
                    flag = true;
                    m_range.erase(std::remove_if(m_range.begin(), m_range.end(), [&](auto && t){return t == r;}), m_range.end());
                    delete r;
                    break;
                }
            }
            if (flag) continue; else break;
        }
        for (auto && range : _fixed_ranges) {
            m_range.push_back(make_range(range[0], range[1], range[0], range[1]));
        }
    }

    void merge_ranges() {
        LOG << "# Merge ranges..." << std::endl;
        while (true) {
            int flag = 0;
            for (auto && r1 : m_range) {
                for (auto && r2 : m_range) {
                    auto &range1 = *r1;
                    auto &range2 = *r2;
                    if (r1 != r2) {
                        if (range1[1] < range1[2]) {
                            if (range1[1] + 1 == range2[0] && range2[3] + 1 == range1[2]) {
                                if (range2[1] < range2[2]) {
                                    m_range.push_back(make_range(range1[0], range2[1], range2[2], range1[3]));
                                } else {
                                    m_range.push_back(make_range(range1[0], range1[3], range1[0], range1[3]));
                                }
                                delete r1;
                                delete r2;
                                m_range.erase(std::remove_if(m_range.begin(), m_range.end(), [&](auto && t){return t == r1 || t == r2;}), m_range.end());
                                flag++;
                                break;
                            }
                        }
                    }
                }
                if (flag > 0) {
                    break;
                }
            }
            if (flag == 0) {
                break;
            }
        }
    }

// MC related methods

    virtual void init_run() {
        LOG << "# Set pseudo-knots" << std::endl;
        set_pseudo_knots();
    }

    virtual void mc_select() {
        int len = m_range.size();
        _moved_residues = *(m_range[int(rand() * len)]);
    }

    virtual bool is_selected(const int &i) const {
        auto &p = _moved_residues;
        return (i >= p[0] && i <= p[1]) || (i >= p[2] && i <= p[3]);
    }

    virtual Vec rotating_center() const {
        int beg = _moved_residues[0];
        int end = _moved_residues[3];
        auto &r1 = mc_t::_pred_chain[beg];
        auto &r2 = mc_t::_pred_chain[end];
        Vec origin = Vec::Zero(3);
        int n_atom = 0;
        for (auto && atom : r1) {
            for (int i = 0; i < 3; i++) origin[i] += atom[i];
            n_atom++;
        }
        for (auto && atom : r2) {
            for (int i = 0; i < 3; i++) origin[i] += atom[i];
            n_atom++;
        }
        for (int i = 0; i < 3; i++) origin[i] /= n_atom;
        return origin;
    }

};

template<typename T>
void chain_refine(Chain &chain, loop *l, const fixed_ranges_t &fixed_ranges = {}, std::string traj = "") {
    Par par;
    par._orig_pars = {"nsp", ""};
    std::string seq, ss;
    seq_read_tree(seq, l);
    ss_read_tree(ss, l);
    par._pars["traj"].push_front(traj);
    par._pars["seq"].push_front(seq);
    par._pars["ss"].push_front(ss);
    DHMC<T> mc(par);
    mc._pred_chain = chain;
    mc._fixed_ranges = fixed_ranges;
    mc.run();
    chain = mc._pred_chain;
}

} // namespace nuc3d
} // namespace jian

