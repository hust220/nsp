#pragma once

#include <Eigen/Dense>
#include <map>
#include <deque>
#include <vector>
#include <numeric>
#include "QHMC.hpp"
#include "Module.hpp"
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../mcsm.hpp"
#include "../utils/Factory.hpp"
#include "../utils/Par.hpp"
#include "../utils/Env.hpp"

BEGIN_JN
namespace qhmc {

using fac_t = Factory<Module::cons_t>;

class QHMC : public MCSM {
public:
    using res_t = struct {char seq; char ss; int num;};
    using res_list_t = std::deque<res_t>;
    using indices_t = std::deque<int>;
    using related_residues_t = std::deque<indices_t *>;

    Tree m_tree;
    std::deque<Module *> m_modules;
    int m_selected_index;
    related_residues_t m_related_residues,
                       m_unrelated_residues,
                       m_all_indices;
    std::array<bool, 4> m_directions; // true: parallel, false: antiparallel
    std::array<int, 4> m_arrangement;
    std::map<indices_t *, Chain> m_fixed_ranges;

    QHMC() = default;

    void init(const Par &par) {
		MCSM::init(par);
    }

    ~QHMC() {
        for (auto && i : m_modules) {
            delete i;
        }
        for (auto && i : m_all_indices) {
            delete i;
        }
    }

    void print_modules() {
        for (auto && module : m_modules) {
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

    void set_modules() {
        int len = _seq.size();
        m_modules.push_back(fac_t::create("head_hairpin", m_tree.front().front(), m_tree.back().back(), len));
        int i = 0;
        for (; i + 1 < m_tree.size(); i++) {
            m_modules.push_back(fac_t::create("helix", m_tree[i].front(), m_tree[i].back(), len));
            m_modules.push_back(fac_t::create("loop", m_tree[i].back(), m_tree[i+1].front(), len));
        }
        m_modules.push_back(fac_t::create("helix", m_tree[i].front(), m_tree[i].back(), len));
        m_modules.push_back(fac_t::create("tail_hairpin", m_tree.front().front(), m_tree.back().back(), len));
    }

    void set_res_list(res_list_t &res_list) {
        int len = _seq.size();
        for (int i = 0; i < len; i++) {
            res_list.push_back({_seq[i], _ss[i], i});
        }
    }

    void ss_to_tree() {
        Tuples tuples;
        tuples_from_ss(tuples);
        tuples_to_tree(tuples);
        print_tree();
    }

    void build_initial_scaffold() {
        LOG << "## Compute maximum length" << std::endl;
        int len = std::accumulate(m_modules.begin(), m_modules.end(), 0, [](int n, auto &&m){
            return n + m->d_max_len;
        });
        LOG << "## Build helix" << std::endl;
        Chain &&c = build_helix(len);
        //mol_write(c, "bb.pdb");
        LOG << "## Shrink to fit" << std::endl;
        shrink_to_fit(c);
        //mol_write(_pred_chain, "cc.pdb");
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
        S file_name = Env::lib() + "/RNA/pars/nuc3d/quadruple/quadruple-helix-" + JN_STR(n) + ".pdb";
        Chain chain;
        chain_read_model(chain, file_name);
        //return cg_t::chain(chain);
        return chain;
    }

    void set_coords_residue(Mat &c1, int m, const Residue &r) {
        static int l = m_cg->res_size();
        Residue res = m_cg->to_cg(r);
        for (int i = 0; i < l; i++) {
            for (int j = 0; j < 3; j++) {
                c1(m * l + i, j) = res[i][j];
            }
        }
    }

    Chain connect_quadruple_helix(Chain &c1, Chain &c2) {
        int l = m_cg->res_size();
        int len1 = c1.size()/4, len2 = c2.size()/4;
        int len = len1 + len2 - 1;
        Mat m1(4*l, 3), m2(4*l, 3);
        Chain c;
        int i, j;

        LOG << "#### Set m1 and m2" << std::endl;
        for (i = 0; i < 4; i++) {
            set_coords_residue(m1, i, c1[len1+i*len1-1]);
            set_coords_residue(m2, i, c2[i*len2]);
        }

        LOG << "#### Supperposition." << std::endl;
        auto sp = geom::suppos(m1, m2);
        INIT_SUPPOS(sp);
		for (auto && res : c1) for (auto && atom : res) { APPLY_SUPPOS(atom, sp); }

        LOG << "#### Set coordinates." << std::endl;
        for (i = 0; i < 4; i++) {
            for (j = 0; j < len1; j++) c.push_back(c1[i*len1+j]);
            for (j = 1; j < len2; j++) c.push_back(c2[i*len2+j]);
        }

        return c;
    }

    void shrink_to_fit(const Chain &c) {
        int i, j, k, l, a, b;

        int len = c.size()/4;
        LOG << "len: " << len << std::endl;
        int n = 0;
        _pred_chain.resize(_seq.size());
        for (i = 0; i < m_modules.size(); i++) {
            Mati &m = m_modules[i]->d_indices;
            l = m.rows();
            LOG << "module: " << i << std::endl;
            LOG << m << std::endl;
            for (j = 0; j < l; j++) {
                for (k = 0; k < 4; k++) {
                    if (m(j, k) != -1) {
                        //a = (k % 2 == 0 ? n+j+k*len : len-1-n-j+k*len);
                        b = std::distance(m_arrangement.begin(), std::find(m_arrangement.begin(), m_arrangement.end(), k));
                        a = b*len+n+j;
                        _pred_chain[m(j, k)] = c[a];
                    }
                }
            }
            n += l;
        }
    }

    void tuples_from_ss(Tuples &tuples) {
        std::deque<int> dq;
        auto & ss = _ss;
        int i, j, n, l, size;

        // Set dq
        i = 0;
        for (auto && c : ss) {
            if (c == 'G') {
                dq.push_back(i);
            } else if (c == 'g') {
                dq.push_back(i);
            } else {
                // ...
            }
            i++;
        }

        size = dq.size();
        l = size / 4;

        // Set directions
        m_directions[0] = true;
        for (i = 1; i < 4; i++) {
            m_directions[i] = (ss[dq[l*i]] == 'G' ? true : false);
        }

        // Set tuples
        Tuple t;
        for (i = 0; i < l; i++) {
            for (j = 0; j < 4; j++) {
                n = (m_directions[j] ? l*j+i : l*j+l-i-1);
                t[j] = dq[n];
            }
            tuples.push_back(t);
        }
    }

    template<typename T, typename U>
    bool adjacent(T &&t1, U &&t2) {
        return abs(t1[0] - t2[0]) == 1 &&
               abs(t1[1] - t2[1]) == 1 &&
               abs(t1[2] - t2[2]) == 1 &&
               abs(t1[3] - t2[3]) == 1;
    }

    void tuples_to_tree(Tuples tuples) {
        Tuples dq;
        dq.push_back(tuples[0]);
        for (int i = 1; i < tuples.size(); i++) {
            if (!(adjacent(tuples[i-1], tuples[i]))) {
                m_tree.push_back(std::move(dq));
            }
            dq.push_back(tuples[i]);
        }
        m_tree.push_back(std::move(dq));
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
        for (auto && helix : m_tree) {
            print_helix(helix);
        }
    }

    // mc-related functions

    void set_unrelated_residues() {
        int len = _seq.size();
        related_residues_t &r = m_related_residues;
        m_unrelated_residues.resize(len);
        for (int i = 0; i < len; i++) {
            m_unrelated_residues[i] = new std::deque<int>();
            for (int j = 0; j < len; j++) {
                if (std::none_of(r[i]->begin(), r[i]->end(), [&j](auto && n){return n == j;})) {
                    m_unrelated_residues[i]->push_back(j);
                }
            }
        }
    }

    void print_related_residues(const related_residues_t &r) {
        int len = _seq.size();
        for (int i = 0; i < len; i++) {
            LOG << i << ' ';
            for (auto && j : *(r[i])) {
                LOG << j << ' ';
            }
            LOG << std::endl;
        }
    }

    void set_related_residues() {
        m_related_residues.resize(_seq.size());
        for (auto && module : m_modules) {
            if (module->type() != "helix") {
                for (auto && frag : module->d_frags) {
                    for (auto && i : frag) {
                        indices_t * p = new indices_t;
                        m_all_indices.push_back(p);
                        m_related_residues[i] = p;
                        m_related_residues[i]->push_back(i);
                    }
                }
            }
        }
        for (auto && module : m_modules) {
            if (module->type() == "helix") {
                indices_t * p = new indices_t;
                m_all_indices.push_back(p);
                for (auto && frag : module->d_frags) {
                    for (auto && i : frag) {
                        m_related_residues[i] = p;
                        p->push_back(i);
                    }
                }
            }
        }
    }

    void set_related_and_unrelated_residues() {
        set_related_residues();
        set_unrelated_residues();
        LOG << "## Print related residues" << std::endl;
        print_related_residues(m_related_residues);
        LOG << "## Print unrelated residues" << std::endl;
        print_related_residues(m_unrelated_residues);
    }

    virtual void before_run() {
        LOG << "# Convert 2D structure to tree." << std::endl;
        ss_to_tree();

        LOG << "# Set modules." << std::endl;
        set_modules();

        LOG << "# Print modules." << std::endl;
        print_modules();

        LOG << "# Build initial scaffold." << std::endl;
        build_initial_scaffold();

        LOG << "# Set related and unrelated residues..." << std::endl;
        set_related_and_unrelated_residues();
    }

    virtual void mc_select() {
        int len = _seq.size();
        m_selected_index = int(rand() * len);
    }


    virtual bool is_selected(const int &i) const {
        auto &v = *(m_related_residues[m_selected_index]);
        return std::find(v.begin(), v.end(), i) != v.end();
    }

    virtual Vec rotating_center() const {
        Vec vec = Vec::Zero(3);
        indices_t &v = *(m_related_residues[m_selected_index]);
        double n = 0;
        for (auto && atom : _pred_chain[v[int(rand()*v.size())]]) {
            for (int j = 0; j < 3; j++) {
                vec[j] += atom[j];
            }
            n++;
        }
        for (int j = 0; j < 3; j++) {
            vec[j] /= n;
        }
        return vec;
    }

    void save_helix() {}

    virtual void save_fixed_ranges() {
        for (auto && indices : m_all_indices) {
            if (indices->size() > 1) {
                Chain c;
                for (auto && i : *indices) {
                    c.push_back(_pred_chain[i]);
                }
                m_fixed_ranges[indices] = c;
            }
        }
    }

    void restore_helix(indices_t * indices) {
        int l = indices->size();
        int s = m_cg->res_size();
        Mat m1(l*s, 3), m2(l*s, 3);
        int n;
        Chain & c = m_fixed_ranges[indices];

        // Set m1 and m2
        n = 0;
        for (auto && i : *indices) {
            set_coords_residue(m1, n, c[n]);
            set_coords_residue(m2, n, _pred_chain[i]);
            n++;
        }

        // Superposition
		geom::Superposition<double> sp(m1, m2);
		for (auto && res : c) for (auto && atom : res) { sp.apply(atom); }

        // Set _pred_chain
        n = 0;
        for (auto && i : *indices) {
            _pred_chain[i] = c[n];
            n++;
        }
    }

    virtual void restore_fixed_ranges() {
        for (auto && indices : m_all_indices) {
            if (indices->size() > 1) {
                restore_helix(indices);
            }
        }
    }

    virtual void read_ss() {
        auto die = [&](){
            throw "jian::QHMC::read_ss error! Illegal secondary structure!";
        };

        auto check_arrangement = [&](auto && s){
            std::set<char> set;
            for (auto && c : s) set.insert(c);
            if (set != std::set<char>{'1', '2', '3', '4'}) die();
        };

        auto check_ss = [&](auto && s) {
            if (!std::regex_match(s, std::regex("^[Gg.-]+$"))) die();
        };

        S ss = _par->get("ss");
        tokenize_v v;
        jian::tokenize(ss, v, ": ");
        std::cout << "v.size() " << v.size() << ' ' << ss << std::endl;
        if (v.size() == 1) {
            check_ss(v[0]);
            _ss = v[0];
        } else if (v.size() == 2) {
            check_arrangement(v[0]);
            for (int i = 0; i < 4; i++) m_arrangement[i] = std::stoi(v[0].substr(i, 1))-1;
            check_ss(v[1]);
            _ss = v[1];
        } else {
            die();
        }
    }

};

} // namespace qhmc
END_JN


