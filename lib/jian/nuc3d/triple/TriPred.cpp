#include <Eigen/Dense>
#include <map>
#include <deque>
#include <vector>
#include "TriPred.hpp"
#include "Module.hpp"
#include "CG2AA.hpp"
#include "../JobPredict3D.hpp" 
#include "../../scoring/new_score.hpp"
#include "../../mc.hpp"
#include "../../pdb.hpp"
#include "../../geom.hpp"
#include "../../utils/Factory.hpp"
#include "../../utils/Par.hpp"
#include "../../utils/rand.hpp"
#include "../../utils/Env.hpp"

namespace jian {
namespace triple {

using fac_t = Factory<Module::cons_t>;

class TriPredImpl : public MC, public JobPredict3D {
public:
    using Res = struct {char seq; char ss; int num;};
    using indices_t = std::deque<int>;
    using related_residues_t = std::vector<std::shared_ptr<indices_t>>;

    Tree _tree;
    Mat d_c;
    std::deque<triple::Module *> d_modules;
    int d_mc_selected_index;
    related_residues_t d_mc_related_residues;
    related_residues_t d_mc_unrelated_residues;
    std::deque<Vec> d_mc_moved_residues;
    std::vector<int> d_indices;
    Chain d_chain;

    double _dist_o3_c5 {3.1};

    std::vector<std::vector<std::string>> _coarse_grained_atoms {
        {"C5*", "O3*", "C1*", "C2", "C6"},
        {"C5*", "O3*", "C1*", "C2", "C4"},
        {"C5*", "O3*", "C1*", "C2", "C6"},
        {"C5*", "O3*", "C1*", "C2", "C4"}
    };

    TriPredImpl() {}

    TriPredImpl(Par par) : JobPredict3D(par), d_c(_seq.size()*5, 3) {
        d_indices.resize(_seq.size());
        for (int i = 0; i < _seq.size(); i++) {
            d_indices[i] = type_id(_seq[i]);
        }
    }

    ~TriPredImpl() {
        for (auto && i : d_modules) {
            delete i;
        }
    }

    void predict() {
        std::cout << "# Convert 2D structure to tree." << std::endl;
        ss_to_tree();
        std::cout << "# Set modules." << std::endl;
        set_modules();
        print_modules();
        std::cout << "# Build initial scaffold." << std::endl;
        build_initial_scaffold();
        mc_init();
        mc();
        cg2aa();
    }

    int type_id(const char &c) {
        if (c == 'A') return 0;
        else if (c == 'T') return 1;
        else if (c == 'G') return 2;
        else if (c == 'C') return 3;
    }

    std::shared_ptr<Mat> load_triple_helix(int n) {
        int num_atoms = n * 3 * 5;
        std::shared_ptr<Mat> c = std::make_shared<Mat>(num_atoms, 3);
        std::string file_name = Env::lib() + "/RNA/pars/nuc3d/triple/triple-helix-" + JN_STR(n) + ".pdb";
        Chain &&chain = residues_from_file(file_name);
        for (int i = 0; i < chain.size(); i++) {
            int type = type_id(chain[i].name.back());
            auto &names = _coarse_grained_atoms[type];
            for (auto && atom : chain[i]) {
                auto r = std::find(names.begin(), names.end(), atom.name);
                if (r != names.end()) {
                    int d = std::distance(names.begin(), r);
                    for (int j = 0; j < 3; j++) {
                        (*c)(i * 5 + d, j) = atom[j];
                    }
                }
            }
        }
        return c;
    }

    void set_coords_residue(Mat &c1, int m, const Mat &c2, int n) {
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 3; j++) {
                c1(m * 5 + i, j) = c2(n * 5 + i, j);
            }
        }
    }

    std::shared_ptr<Mat> connect_triple_helix(Mat &c1, Mat &c2) {
        std::cout << "c1: \n" << c1 << std::endl;
        std::cout << "c2: \n" << c2 << std::endl;
        int len1 = c1.rows()/15, len2 = c2.rows()/15;
        int len = len1 + len2 - 1;
        Mat m1(15, 3), m2(15, 3);
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 3; j++) {
                m1(i, j) = c1((len1 - 1) * 5 + i, j);
                m1(5 + i, j) = c1((len1) * 5 + i, j);
                m1(10 + i, j) = c1((3 * len1 - 1) * 5 + i, j);
                m2(i, j) = c2(i, j);
                m2(5 + i, j) = c2((2 * len2 - 1) * 5 + i, j);
                m2(10 + i, j) = c2((2 * len2) * 5 + i, j);
            }
        }
        std::cout << "m1: \n" << m1 << std::endl;
        std::cout << "m2: \n" << m2 << std::endl;
        std::cout << "#### Supperposition." << std::endl;
        auto sp = geom::suppos(m2, m1);
        INIT_SUPPOS(sp);
        APPLY_SUPPOS_m(c2, sp);
        std::shared_ptr<Mat> c = std::make_shared<Mat>(len*15, 3);
        std::cout << "#### Set coordinates." << std::endl;
        for (int i = 0; i < len1; i++) set_coords_residue(*c, i, c1, i);
        for (int i = 0; i < len2 * 2 - 2; i++) set_coords_residue(*c, len1 + i, c2, i + 1);
        for (int i = 0; i < len1 * 2; i++) set_coords_residue(*c, len1 + len2 * 2 - 2 + i, c1, len1 + i);
        for (int i = 0; i < len2 - 1; i++) set_coords_residue(*c, len1 * 3 + len2 * 2 - 2 + i, c2, 2 * len2 + 1 + i);
        return c;
    }

    void mc_write() {
//        static int n = 1;
//        if (n == 1) {
//            file::clear_file(_name + ".tripred.mc.pdb");
//        }
//        append_chain_to_file(_pred_chain, _name + ".tripred.mc.pdb", n);
        std::cout << _mc_step + 1 << ": "<< _mc_en << "(energy) " << _mc_tempr << "(temprature) " << _mc_local_succ_rate << "(success rate)" << std::endl;
//        n++;
    }

    Vec center(int n) {
        Vec v = Vec::Zero(3);
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 3; j++) {
                v[j] += d_c(n * 5 + i, j);
            }
        }
        return v / 5.0;
    }

    double mc_partial_energy() {
        double e, d;
        indices_t &v1 = *(d_mc_related_residues[d_mc_selected_index]);   
        indices_t &v2 = *(d_mc_unrelated_residues[d_mc_selected_index]);   
        Eigen::MatrixXd m(10, 3);
        for (auto && i : v1) {
            for (auto && j : v2) {
                d = geom::distance(center(i), center(j));
                if (d < 20) {
                    set_coords_residue(m, 0, d_c, i);
                    set_coords_residue(m, 1, d_c, j);
                    e += scoring::new_score(m, d_indices[i] * 4 + d_indices[j]);
                }
                if ((i - j == 1 || j - i == 1) && d > 6) {
                    e += 9999;
                }
            }
        }
        return e;
    }

    void mc_select() {
        int len = _seq.size();
        d_mc_selected_index = int(rand() * len);
        mc_backup();
    }

    Vec center_selected_residues() {
        Vec vec = Vec::Zero(3);
        int len = d_mc_moved_residues.size();
        for (int i = 0; i < len; i++) {
            for (int j = 0; j < 3; j++) {
                vec[j] += d_c(d_mc_related_residues[d_mc_selected_index]->at(i), j);
            }
        }
        for (int j = 0; j < 3; j++) {
            vec[j] /= len;
        }
        return vec;
    }

    void mc_sample() {
        int len = d_mc_moved_residues.size();
        if (rand() < 0.5) {
            // translate
            int index = int(rand() * 3);
            double dist = (rand() - 0.5) * 2;
            for (int i = 0; i < len; i++) {
//                std::cout << i << ' ' << len << std::endl;
                d_c(d_mc_related_residues[d_mc_selected_index]->at(i), index) += dist;
            }
        } else {
            // rotate
            int index = int(rand() * 3);
            double dih = (rand() - 0.5) * PI / 6;
            Mat &&rot = geom::rot_mat(index, dih);
            Vec &&origin = center_selected_residues();
            for (int i = 0; i < len; i++) {
                geom::rotate(d_c.row(d_mc_related_residues[d_mc_selected_index]->at(i)), origin, rot);
            }
        }
    }

    void mc_backup() {
        d_mc_moved_residues.clear();
        int len = d_mc_related_residues[d_mc_selected_index]->size();
        for (int i = 0; i < len; i++) {
            d_mc_moved_residues.push_back(d_c.row(d_mc_related_residues[d_mc_selected_index]->at(i)));
        }
    }

    void mc_back() {
        int len = d_mc_moved_residues.size();
        for (int i = 0; i < len; i++) {
            d_c.row(d_mc_related_residues[d_mc_selected_index]->at(i)) = d_mc_moved_residues[i];
        }
    }

    void mc_init() {
        d_mc_related_residues.resize(_seq.size());
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
        print_related_residues(d_mc_related_residues);
        print_related_residues(d_mc_unrelated_residues);
    }

    void set_unrelated_residues() {
        related_residues_t &r = d_mc_related_residues;
        d_mc_unrelated_residues.resize(_seq.size());
        for (int i = 0; i < _seq.size(); i++) {
            d_mc_unrelated_residues[i] = std::make_shared<std::deque<int>>();
            for (int j = 0; j < _seq.size(); j++) {
                if (std::none_of(r[i]->begin(), r[i]->end(), [&j](auto && n){return n == j;})) {
                    d_mc_unrelated_residues[i]->push_back(j);
                }
            }
        }
    }

    void print_related_residues(const related_residues_t &r) {
        for (int i = 0; i < _seq.size(); i++) {
            std::cout << i << ' ';
            for (auto && j : *(r[i])) {
                std::cout << j << ' ';
            }
            std::cout << std::endl;
        }
    }

    void cg2aa() {
        std::vector<std::string> v(_seq.size());
        for (int i = 0; i < v.size(); i++) {
            v[i] = std::string("D") + _seq[i];
        }
        d_chain = cg_to_aa(d_c, v);
    }

    void print_modules() {
        for (auto && module : d_modules) {
            std::cout << module << ' ' << module->d_max_len;
            for (auto && frag : module->d_frags) {
                std::cout << ' ';
                for (auto && i : frag) {
                    std::cout << i << '-';
                }
            }
            std::cout << std::endl;
        }
    }

    void set_modules() {
        d_modules.push_back(fac_t::create("head_hairpin", _tree.front().front(), Tuple{0, int(_seq.size()), 0}));
        int i = 0;
        for (; i + 1 < _tree.size(); i++) {
            d_modules.push_back(fac_t::create("helix", _tree[i].front(), _tree[i].back()));
            d_modules.push_back(fac_t::create("loop", _tree[i].back(), _tree[i+1].front()));
        }
        d_modules.push_back(fac_t::create("helix", _tree[i].front(), _tree[i].back()));
        d_modules.push_back(fac_t::create("tail_hairpin", _tree.back().back(), Tuple{0, int(_seq.size()), 0}));
    }

    void ss_to_tree() {
        int len = _seq.size();
        std::deque<Res> res_list;
        for (int i = 0; i < len; i++) {
            res_list.push_back({_seq[i], _ss[i], i});
        }
        Tuples &&tuples = get_tuples(res_list);
        print_helix(tuples);
        tuples_to_tree(tuples);
        print_tree();
    }

    void build_initial_scaffold() {
        std::cout << "## Compute maximum length." << std::endl;
        int len = std::accumulate(d_modules.begin(), d_modules.end(), 0, [](int n, auto &&m){
            return n + m->d_max_len;
        });
        std::cout << "## Build helix." << std::endl;
        std::shared_ptr<Mat> c = build_helix(len);
        std::cout << "## Shrink to fit." << std::endl;
        shrink_to_fit(*c);
        std::cout << d_c << std::endl;
    }

    std::shared_ptr<Mat> build_helix(int len) {
        std::shared_ptr<Mat> c, c_;
        if (len <= 2) {
            std::cout << "### Load triple helix." << std::endl;
            c = load_triple_helix(len);
        } else {
            std::cout << "### Load triple helix." << std::endl;
            c = load_triple_helix(2);
            for (int i = 2; i < len; i++) {
                std::cout << "### Load triple helix." << std::endl;
                c_ = load_triple_helix(2);
                std::cout << "### Connect triple helix." << std::endl;
                c = connect_triple_helix(*c, *c_);
            }
        }
        return c;
    }

    void shrink_to_fit(const Mat &c) {
        std::cout << d_c.rows() << std::endl;
        int len = c.rows()/15;
        std::cout << "len: " << len << std::endl;
        int n = 0;
        for (int i = 0; i < d_modules.size(); i++) {
            Mat &m = *(d_modules[i]->d_indices);
            int l = m.rows();
            std::cout << m << std::endl;
            for (int j = 0; j < l; j++) {
                std::cout << i << ' ' << j << std::endl;
                if (m(j, 0) != -1) {
                    set_coords_residue(d_c, m(j, 0), c, n + j);
                }
                if (m(j, 1) != -1) {
                    set_coords_residue(d_c, m(j, 1), c, 2 * len - 1 - n - j);
                }
                if (m(j, 2) != -1) {
                    set_coords_residue(d_c, m(j, 2), c, 2 * len + n + j);
                }
            }
            n += l;
        }
    }

    template<typename T>
    Tuples get_tuples(T &&res_list) {
        std::deque<std::array<int, 3>> tuples;
        std::vector<int> v(res_list.size());
        char b;
        int flag = 0;
        for (int i = 0; i < res_list.size(); i++) {
            b = res_list[i].ss;
            if (b == '1' || b == '2' || b == '3') {
                v[flag] = i;
                flag++;
            }
        }
        int len = (flag+1)/3;
        Tuple t;
        for (int i = 0; i < len; i++) {
            t = {v[i], v[2*len-1-i], v[2*len+i]};
            std::sort(t.begin(), t.end(), [&res_list](int a, int b){return res_list[a].ss < res_list[b].ss;});
            tuples.push_back(std::move(t));
        }
        return tuples;
    }

    template<typename T, typename U>
    bool adjacent(T &&t1, U &&t2) {
        return abs(t1[0] - t2[0]) == 1 &&
               abs(t1[1] - t2[1]) == 1 &&
               abs(t1[2] - t2[2]) == 1;
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
        std::cout << tuple[0] << ' ' << tuple[1] << ' ' << tuple[2] << std::endl;
    }

    void print_helix(const Tuples &helix) {
        std::cout << "Helix:" << std::endl;
        for (auto && tuple : helix) {
            print_tuple(tuple);
        }
    }

    void print_tree() {
        std::cout << "Tree: " << std::endl;
        for (auto && helix : _tree) {
            print_helix(helix);
        }
    }

};

} // namespace triple

TriPred::TriPred(const Par &par) : _impl(std::make_shared<triple::TriPredImpl>(par)) {}

void TriPred::predict() {
    _impl->predict();
}

} // namespace jian


