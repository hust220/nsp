#pragma once

#include <Eigen/Dense>
#include <map>
#include <deque>
#include <vector>
#include "QHMC.hpp"
#include "Module.hpp"
#include "../mc/MC1p.hpp"
#include "../mc/MCpsb.hpp"
#include "../../pdb.hpp"
#include "../../geom.hpp"
#include "../../utils/Factory.hpp"
#include "../../utils/Par.hpp"
#include "../../utils/Env.hpp"

namespace jian {
namespace nuc3d {
namespace quadruple {

using fac_t = Factory<Module::cons_t>;

template<typename MC_T>
class QHMC : public MC_T {
public:
    using mc_t = MC_T;

    using Res = struct {char seq; char ss; int num;};
    using indices_t = std::deque<int>;
    using related_residues_t = std::vector<std::shared_ptr<indices_t>>;

    Tree _tree;
    std::deque<Module *> d_modules;
    int d_mc_selected_index;
    related_residues_t d_mc_related_residues;
    related_residues_t d_mc_unrelated_residues;

    QHMC(Par par) : mc_t(par) {}

    ~QHMC() {
        for (auto && i : d_modules) {
            delete i;
        }
    }

    void set_modules() {
        int len = mc_t::_seq.size();
        d_modules.push_back(fac_t::create("head_hairpin", _tree.front().front(), Tuple{0, len, 0, 0}));
        int i = 0;
        for (; i + 1 < _tree.size(); i++) {
            d_modules.push_back(fac_t::create("helix", _tree[i].front(), _tree[i].back()));
            d_modules.push_back(fac_t::create("loop", _tree[i].back(), _tree[i+1].front()));
        }
        d_modules.push_back(fac_t::create("helix", _tree[i].front(), _tree[i].back()));
        d_modules.push_back(fac_t::create("tail_hairpin", _tree.back().back(), Tuple{0, len, 0, 0}));
        for (auto && module : d_modules) {
            LOG << module << ' ' << module->d_max_len;
            for (auto && frag : module->d_frags) {
                LOG << ' ';
                for (auto && i : frag) {
                    LOG << i << '-';
                }
            }
            LOG << std::endl;
        }
    }

    void ss_to_tree() {
        int len = mc_t::_seq.size();
        std::deque<Res> res_list;
        for (int i = 0; i < len; i++) {
            res_list.push_back({mc_t::_seq[i], mc_t::_ss[i], i});
        }
        Tuples &&tuples = get_tuples(res_list);
        print_helix(tuples);
        tuples_to_tree(tuples);
        print_tree();
    }

    void build_initial_scaffold() {
        LOG << "## Compute maximum length." << std::endl;
        int len = std::accumulate(d_modules.begin(), d_modules.end(), 0, [](int n, auto &&m){
            return n + m->d_max_len;
        });
        LOG << "## Build helix." << std::endl;
        Chain &&c = build_helix(len);
        mol_write(c, "bb.pdb");
        LOG << "## Shrink to fit." << std::endl;
        shrink_to_fit(c);
    }

    Chain build_helix(int len) {
        Chain c, c_;
        if (len <= 2) {
            LOG << "### Load quadruple helix." << std::endl;
            c = load_quadruple_helix(len);
        } else {
            LOG << "### Load quadruple helix." << std::endl;
            c = load_quadruple_helix(2);
            for (int i = 2; i < len; i++) {
                LOG << "### Load quadruple helix." << std::endl;
                c_ = load_quadruple_helix(2);
                LOG << "### Connect quadruple helix." << std::endl;
                c = connect_quadruple_helix(c, c_);
            }
        }
        return c;
    }

    Chain load_quadruple_helix(int n) {
        std::string file_name = Env::lib() + "/RNA/pars/nuc3d/quadruple/quadruple-helix-" + JN_STR(n) + ".pdb";
        Chain chain;
        chain_read_model(chain, file_name);
        return mc_t::cg_t::chain(chain);
    }

    void set_coords_residue(Mat &c1, int m, const Residue &r) {
        static int l = mc_t::cg_t::size_res;
        for (int i = 0; i < l; i++) {
            for (int j = 0; j < 3; j++) {
                c1(m * l + i, j) = r[i][j];
            }
        }
    }

    Chain connect_quadruple_helix(Chain &c1, Chain &c2) {
        int l = mc_t::cg_t::size_res;
        int len1 = c1.size()/4, len2 = c2.size()/4;
        int len = len1 + len2 - 1;
        Mat m1(4*l, 3), m2(4*l, 3);
        set_coords_residue(m1, 0, c1[len1-1]);
        set_coords_residue(m1, 1, c1[len1]);
        set_coords_residue(m1, 2, c1[3*len1-1]);
        set_coords_residue(m1, 3, c1[3*len1]);
        set_coords_residue(m2, 0, c2[0]);
        set_coords_residue(m2, 1, c2[2*len2-1]);
        set_coords_residue(m2, 2, c2[2*len2]);
        set_coords_residue(m2, 3, c2[4*len2-1]);
        LOG << "#### Supperposition." << std::endl;
        auto sp = geom::suppos(m1, m2);
        INIT_SUPPOS(sp);
        for (auto && res : c1) for (auto && atom : res) APPLY_SUPPOS(atom, sp);
        Chain c;
        LOG << "#### Set coordinates." << std::endl;
        for (int i = 0; i < len1; i++)     c.push_back(c1[i]);
        for (int i = 0; i < 2*len2-2; i++) c.push_back(c2[1+i]);
        for (int i = 0; i < 2*len1; i++)   c.push_back(c1[len1+i]);
        for (int i = 0; i < 2*len2-2; i++) c.push_back(c2[2*len2+1+i]);
        for (int i = 0; i < len1; i++)     c.push_back(c1[3*len1+i]);
        return c;
    }

    void shrink_to_fit(const Chain &c) {
        int len = c.size()/4;
        LOG << "len: " << len << std::endl;
        int n = 0;
        for (int i = 0; i < d_modules.size(); i++) {
            Mat &m = *(d_modules[i]->d_indices);
            int l = m.rows();
            LOG << m << std::endl;
            for (int j = 0; j < l; j++) {
                LOG << i << ' ' << j << std::endl;
                if (m(j, 0) != -1) {
                    mc_t::_pred_chain[m(j, 0)] = c[n+j];
                }
                if (m(j, 1) != -1) {
                    mc_t::_pred_chain[m(j, 1)] = c[2*len-1-n-j];
                }
                if (m(j, 2) != -1) {
                    mc_t::_pred_chain[m(j, 2)] = c[2*len+n+j];
                }
                if (m(j, 3) != -1) {
                    mc_t::_pred_chain[m(j, 3)] = c[4*len-1-n-j];
                }
            }
            n += l;
        }
    }

    template<typename T>
    Tuples get_tuples(T &&res_list) {
        Tuples tuples;
        std::vector<int> v(res_list.size());
        char b;
        int flag = 0;
        for (int i = 0; i < res_list.size(); i++) {
            b = res_list[i].ss;
            if (b == '1') {
                v[flag] = i;
                flag++;
            }
        }
        int len = (flag+1)/4;
        Tuple t;
        for (int i = 0; i < len; i++) {
            t = {v[i], v[2*len-1-i], v[2*len+i], v[4*len-1-i]};
            std::sort(t.begin(), t.end(), [&res_list](int a, int b){return res_list[a].ss < res_list[b].ss;});
            tuples.push_back(std::move(t));
        }
        return tuples;
    }

    template<typename T, typename U>
    bool adjacent(T &&t1, U &&t2) {
        return abs(t1[0] - t2[0]) == 1 &&
               abs(t1[1] - t2[1]) == 1 &&
               abs(t1[2] - t2[2]) == 1 &&
               abs(t1[3] - t2[3]) == 1;
    }

    template<typename T>
    void tuples_to_tree(T &&tuples) {
        std::sort(tuples.begin(), tuples.end(),
            [](auto &&tuple1, auto &&tuple2) {
                return *(std::min_element(tuple1.begin(), tuple1.end())) < 
                       *(std::min_element(tuple2.begin(), tuple2.end()));
            }
        );
        Tuples dq;
        dq.push_back(tuples[0]);
        for (int i = 1; i < tuples.size(); i++) {
            if (!(adjacent(tuples[i-1], tuples[i]))) {
                _tree.push_back(std::move(dq));
            }
            dq.push_back(tuples[i]);
        }
        _tree.push_back(std::move(dq));
    }

    void print_tuple(const Tuple &tuple) {
        LOG << tuple[0] << ' ' << tuple[1] << ' ' << tuple[2] << ' ' << tuple[3] << std::endl;
    }

    void print_helix(const Tuples &helix) {
        LOG << "Helix:" << std::endl;
        for (auto && tuple : helix) {
            print_tuple(tuple);
        }
    }

    void print_tree() {
        LOG << "Tree: " << std::endl;
        for (auto && helix : _tree) {
            print_helix(helix);
        }
    }

    // mc-related functions

    void set_unrelated_residues() {
        int len = mc_t::_seq.size();
        related_residues_t &r = d_mc_related_residues;
        d_mc_unrelated_residues.resize(len);
        for (int i = 0; i < len; i++) {
            d_mc_unrelated_residues[i] = std::make_shared<std::deque<int>>();
            for (int j = 0; j < len; j++) {
                if (std::none_of(r[i]->begin(), r[i]->end(), [&j](auto && n){return n == j;})) {
                    d_mc_unrelated_residues[i]->push_back(j);
                }
            }
        }
    }

    void print_related_residues(const related_residues_t &r) {
        int len = mc_t::_seq.size();
        for (int i = 0; i < len; i++) {
            LOG << i << ' ';
            for (auto && j : *(r[i])) {
                LOG << j << ' ';
            }
            LOG << std::endl;
        }
    }

    void mc_init() {
        d_mc_related_residues.resize(mc_t::_seq.size());
        for (auto && module : d_modules) {
            if (module->type() != "helix") {
                for (auto && frag : module->d_frags) {
                    for (auto && i : frag) {
                        d_mc_related_residues[i] = std::make_shared<std::deque<int>>();
                        d_mc_related_residues[i]->push_back(i);
                    }
                }
            }
        }
        for (auto && module : d_modules) {
            if (module->type() == "helix") {
                std::shared_ptr<std::deque<int>> p = std::make_shared<std::deque<int>>();
                for (auto && frag : module->d_frags) {
                    for (auto && i : frag) {
                        d_mc_related_residues[i] = p;
                        p->push_back(i);
                    }
                }
            }
        }
        set_unrelated_residues();
        LOG << "## Print related residues" << std::endl;
        print_related_residues(d_mc_related_residues);
        LOG << "## Print unrelated residues" << std::endl;
        print_related_residues(d_mc_unrelated_residues);
    }

    virtual void init_run() {
        LOG << "# Convert 2D structure to tree." << std::endl;
        ss_to_tree();

        LOG << "# Set modules." << std::endl;
        set_modules();

        LOG << "# Build initial scaffold." << std::endl;
        build_initial_scaffold();

        LOG << "# MC initialization..." << std::endl;
        mc_init();
    }

    virtual void mc_select() {
        int len = mc_t::_seq.size();
        d_mc_selected_index = int(rand() * len);
    }


    virtual bool is_selected(const int &i) const {
        auto &v = *(d_mc_related_residues[d_mc_selected_index]);
        return std::find(v.begin(), v.end(), i) != v.end();
    }

    virtual Vec rotating_center() const {
        Vec vec = Vec::Zero(3);
        indices_t &v = *(d_mc_related_residues[d_mc_selected_index]);
        double n = 0;
//        for (int i = 0; i < v.size(); i++) {
            for (auto && atom : mc_t::_pred_chain[v[int(rand()*v.size())]]) {
                for (int j = 0; j < 3; j++) {
                    vec[j] += atom[j];
                }
                n++;
            }
//        }
        for (int j = 0; j < 3; j++) {
            vec[j] /= n;
        }
        return vec;
    }

};

} // namespace quadruple
} // namespace nuc3d
} // namespace jian


