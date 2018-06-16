#include "g4dna_tsp.hpp"

namespace jian {

namespace qhmc {

    void QHMC::init(const Par &par) {
        MCSM::init(par);
    }

    QHMC::~QHMC() {
        for (auto && i : m_modules) {
            delete i;
        }
        for (auto && i : m_all_indices) {
            delete i;
        }
    }

    void QHMC::print_modules() {
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

    void QHMC::set_modules() {
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

    void QHMC::set_res_list(res_list_t &res_list) {
        int len = _seq.size();
        for (int i = 0; i < len; i++) {
            res_list.push_back({_seq[i], _ss[i], i});
        }
    }

    void QHMC::ss_to_tree() {
        Tuples tuples;
        tuples_from_ss(tuples);
        tuples_to_tree(tuples);
        print_tree();
    }

    void QHMC::build_initial_scaffold() {
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

    Chain QHMC::build_helix(int len) {
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

    Chain QHMC::load_quadruple_helix(int n) {
        S file_name = Env::lib() + "/RNA/pars/nuc3d/quadruple/quadruple-helix-" + JN_STR(n) + ".pdb";
        Chain chain;
        chain_read_model(chain, file_name);
        //return cg_t::chain(chain);
        return chain;
    }

    void QHMC::set_coords_residue(Mat &c1, int m, const Residue &r) {
        static int l = m_cg->res_size();
        Residue res = m_cg->to_cg(r);
        for (int i = 0; i < l; i++) {
            for (int j = 0; j < 3; j++) {
                c1(m * l + i, j) = res[i][j];
            }
        }
    }

    Chain QHMC::connect_quadruple_helix(Chain &c1, Chain &c2) {
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
        for (auto && res : c1) for (auto && atom : res) sp.apply(atom);

        LOG << "#### Set coordinates." << std::endl;
        for (i = 0; i < 4; i++) {
            for (j = 0; j < len1; j++) c.push_back(c1[i*len1+j]);
            for (j = 1; j < len2; j++) c.push_back(c2[i*len2+j]);
        }

        return c;
    }

    void QHMC::shrink_to_fit(const Chain &c) {
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

    void QHMC::tuples_from_ss(Tuples &tuples) {
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

    void QHMC::tuples_to_tree(Tuples tuples) {
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

    void QHMC::print_tuple(const Tuple &tuple) {
        LOG << tuple[0] << ' ' << tuple[1] << ' ' << tuple[2] << ' ' << tuple[3] << std::endl;
    }

    void QHMC::print_helix(const Tuples &helix) {
        LOG << "Helix:" << std::endl;
        for (auto && tuple : helix) {
            print_tuple(tuple);
        }
    }

    void QHMC::print_tree() {
        LOG << "Tree: " << std::endl;
        for (auto && helix : m_tree) {
            print_helix(helix);
        }
    }

    // mc-related functions

    void QHMC::set_unrelated_residues() {
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

    void QHMC::print_related_residues(const related_residues_t &r) {
        int len = _seq.size();
        for (int i = 0; i < len; i++) {
            LOG << i << ' ';
            for (auto && j : *(r[i])) {
                LOG << j << ' ';
            }
            LOG << std::endl;
        }
    }

    void QHMC::set_related_residues() {
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

    void QHMC::set_related_and_unrelated_residues() {
        set_related_residues();
        set_unrelated_residues();
        LOG << "## Print related residues" << std::endl;
        print_related_residues(m_related_residues);
        LOG << "## Print unrelated residues" << std::endl;
        print_related_residues(m_unrelated_residues);
    }

    void QHMC::before_run() {
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

    void QHMC::mc_sample() {
    }

    void QHMC::mc_rollback() {
    }

    void QHMC::mc_backup() {
    }

    void QHMC::mc_select() {
        int len = _seq.size();
        m_selected_index = int(rand() * len);
    }


    bool QHMC::is_selected(const int &i) const {
        auto &v = *(m_related_residues[m_selected_index]);
        return std::find(v.begin(), v.end(), i) != v.end();
    }

    Vec QHMC::rotating_center() const {
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

    void QHMC::save_helix() {}

    void QHMC::save_fixed_ranges() {
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

    void QHMC::restore_helix(indices_t * indices) {
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

    void QHMC::restore_fixed_ranges() {
        for (auto && indices : m_all_indices) {
            if (indices->size() > 1) {
                restore_helix(indices);
            }
        }
    }

    void QHMC::read_ss() {
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

} // namespace qhmc

struct TspG4 {
    Str seq;
    Str ss;
    Chain pred;
    Int n_layers;
};

static Vector<Char> g4_symbols {'A', 'B', 'C', 'D', '.'};

#define JN_DIE(msg) throw to_str("Error in ", __FILE__, ':', __LINE__, "\n", msg)

static int tsp_g4_check_ss(Str ss_) {
    // Ensure that the secondary structure has only the 5 characters: A, U, G, C, .
    Str ss = to_upper_copy(ss_);
    auto & v = g4_symbols;
    if (std::find_if(ss.begin(), ss.end(), [&v](auto &c){ return std::find(v.begin(), v.end(), c)==v.end(); }) != ss.end()) JN_DIE("");
    // Ensure that the number of residues per chain of G4 is the same
    Array<Int, 4> a;
    for (Int i = 0; i < 4; i++) {
        a[i] = std::count(ss.begin(), ss.end(), v[i]);
    }
    if (a[0] != a[1] || a[1] != a[2] || a[2] != a[3]) JN_DIE("");
    return a[0];
}

static void tsp_g4_read_pars(TspG4 &tsp, const Par &par) {
    par.set(tsp.seq, "seq");
    par.set(tsp.ss, "ss");
    if (size(tsp.seq) != size(tsp.ss)) JN_DIE("The length of the sequence and the secondary structure should be the same!");
    tsp.n_layers = tsp_g4_check_ss(tsp.ss);
}

static Chain build_g4(Int n) {
    return Chain{};
}

static Vector<Int> g4_ss_inds(Str ss) {
    auto & symb = g4_symbols;

    Vector<Deque<Int>> inds(symb.size());
    Int n = size(ss);

    Vector<Int> v(n);
    for (Int i = 0; i < n; i++) v[i] = -1;

    for (Int i = 0; i < n; i++) {
        Char c = ss[i];
        auto it = std::find(symb.begin(), symb.begin(), std::toupper(c));
        if (it != symb.end()) {
            Int d = std::distance(symb.begin(), it);
            if (c == std::toupper(c)) inds[d].push_back(i);
            else if (c == std::tolower(c)) inds[d].push_front(i);
        }
    }

    return v;
}

static void tsp_g4_build_scaffold(TspG4 &tsp) {
    Int n = tsp.seq.size();
    Chain g4 = build_g4(tsp.n_layers);
    Vector<Int> v = g4_ss_inds(tsp.ss);
    tsp.pred.resize(n);
    for (Int i = 0; i < n; i++) {
        if (v[i] != -1) tsp.pred[i] = g4[v[i]];
    }
}

static void tsp_g4_opt(TspG4 &tsp) {
}

void tsp_g4(const Par &par) {
    TspG4 tsp;
    tsp_g4_read_pars(tsp, par);
    tsp_g4_build_scaffold(tsp);
    tsp_g4_opt(tsp);
}

}


