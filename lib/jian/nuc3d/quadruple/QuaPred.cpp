#include <Eigen/Dense>
#include <map>
#include <deque>
#include <vector>
#include "QuaPred.hpp"
#include "Module.hpp"
#include "../JobPredict3D.hpp" 
#include "../../mc.hpp"
#include "../../pdb.hpp"
#include "../../geom.hpp"
#include "../../utils/Factory.hpp"
#include "../../utils/Par.hpp"
#include "../../utils/Env.hpp"

namespace jian {
namespace quadruple {

using fac_t = Factory<Module::cons_t>;

class QuaPredImpl : public MC, public JobPredict3D {
public:
    using Res = struct {char seq; char ss; int num;};

    Tree _tree;
    Mat d_c;
    std::deque<quadruple::Module *> d_modules;

    double _dist_o3_c5 {3.1};

    std::vector<std::vector<std::string>> _coarse_grained_atoms {
        {"C5*", "O3*", "C1*", "C2", "C6"},
        {"C5*", "O3*", "C1*", "C2", "C4"},
        {"C5*", "O3*", "C1*", "C2", "C6"},
        {"C5*", "O3*", "C1*", "C2", "C4"}
    };

    QuaPredImpl() {}

    QuaPredImpl(Par par) : JobPredict3D(par), d_c(_seq.size()*5, 3) {}

    ~QuaPredImpl() {
        for (auto && i : d_modules) {
            delete i;
        }
    }

    int type_id(const char &c) {
        if (c == 'A') return 0;
        else if (c == 'T') return 1;
        else if (c == 'G') return 2;
        else if (c == 'C') return 3;
    }

    std::shared_ptr<Mat> load_quadruple_helix(int n) {
        int num_atoms = n * 4 * 5;
        std::shared_ptr<Mat> c = std::make_shared<Mat>(num_atoms, 3);
        std::string file_name = Env::lib() + "/RNA/pars/nuc3d/quadruple/quadruple-helix-" + JN_STR(n) + ".pdb";
        auto &&chain = residues_from_file(file_name);
        for (int i = 0; i < chain.size(); i++) {
            int type = type_id(chain[i].name[1]);
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
//        std::cout << __FUNCTION__ << std::endl;
//        std::cout << c1.rows() << ' ' << c2.rows() << std::endl;
//        std::cout << m << ' ' << n << std::endl;
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 3; j++) {
                c1(m * 5 + i, j) = c2(n * 5 + i, j);
            }
        }
    }

    std::shared_ptr<Mat> connect_quadruple_helix(Mat &c1, Mat &c2) {
        std::cout << "c1: \n" << c1 << std::endl;
        std::cout << "c2: \n" << c2 << std::endl;
        int len1 = c1.rows()/20, len2 = c2.rows()/20;
        int len = len1 + len2 - 1;
        Mat m1(20, 3), m2(20, 3);
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 3; j++) {
                m1(i, j) = c1((len1 - 1) * 5 + i, j);
                m1(5 + i, j) = c1((len1) * 5 + i, j);
                m1(10 + i, j) = c1((3 * len1 - 1) * 5 + i, j);
                m1(15 + i, j) = c1((3 * len1) * 5 + i, j);
                m2(i, j) = c2(i, j);
                m2(5 + i, j) = c2((2 * len2 - 1) * 5 + i, j);
                m2(10 + i, j) = c2((2 * len2) * 5 + i, j);
                m2(15 + i, j) = c2((4 * len2 - 1) * 5 + i, j);
            }
        }
        std::cout << "m1: \n" << m1 << std::endl;
        std::cout << "m2: \n" << m2 << std::endl;
        std::cout << "#### Supperposition." << std::endl;
        auto sp = geom::suppos(m2, m1);
        INIT_SUPPOS(sp);
        APPLY_SUPPOS_m(c2, sp);
        std::shared_ptr<Mat> c = std::make_shared<Mat>(len*20, 3);
        std::cout << "#### Set coordinates." << std::endl;
        int a;
        for (int i = 0; i < len1 * 4; i++) {
            a = int((i/len1+1)/2)*2;
            set_coords_residue(*c, a*(len2-1)+i, c1, i);
        }
        for (int i = 0; i < (len2-1)*4; i++) {
            a = int((i/(len2-1))/2)*2+1;
            set_coords_residue(*c, a*len1+i, c1, a+i);
        }
        return c;
    }

    void mc_write() {}

    double mc_partial_energy() {}

    void mc_select() {}

    void mc_back() {}

    void predict() {
        std::cout << "# Convert 2D structure to tree." << std::endl;
        ss_to_tree();
        std::cout << "# Set modules." << std::endl;
        set_modules();
        std::cout << "# Build initial scaffold." << std::endl;
        build_initial_scaffold();
        return;
        mc();
    }

    void set_modules() {
        d_modules.push_back(fac_t::create("head_hairpin", _tree.front().front(), Tuple{0, int(_seq.size()), 0, 0}));
        int i = 0;
        for (; i + 1 < _tree.size(); i++) {
            d_modules.push_back(fac_t::create("helix", _tree[i].front(), _tree[i].back()));
            d_modules.push_back(fac_t::create("loop", _tree[i].back(), _tree[i+1].front()));
        }
        d_modules.push_back(fac_t::create("helix", _tree[i].front(), _tree[i].back()));
        d_modules.push_back(fac_t::create("tail_hairpin", _tree.back().back(), Tuple{0, int(_seq.size()), 0, 0}));
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
            std::cout << "### Load quadruple helix." << std::endl;
            c = load_quadruple_helix(len);
        } else {
            std::cout << "### Load quadruple helix." << std::endl;
            c = load_quadruple_helix(2);
            for (int i = 2; i < len; i++) {
                std::cout << "### Load quadruple helix." << std::endl;
                c_ = load_quadruple_helix(2);
                std::cout << "### Connect quadruple helix." << std::endl;
                c = connect_quadruple_helix(*c, *c_);
            }
        }
        return c;
    }

    void shrink_to_fit(const Mat &c) {
        std::cout << d_c.rows() << std::endl;
        int len = c.rows()/20;
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
                if (m(j, 3) != -1) {
                    set_coords_residue(d_c, m(j, 3), c, 4 * len - 1 - n - j);
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
        std::cout << tuple[0] << ' ' << tuple[1] << ' ' << tuple[2] << ' ' << tuple[3] << std::endl;
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

} // namespace quadruple

QuaPred::QuaPred(const Par &par) : _impl(std::make_shared<quadruple::QuaPredImpl>(par)) {}

void QuaPred::predict() {
    _impl->predict();
}

} // namespace jian


