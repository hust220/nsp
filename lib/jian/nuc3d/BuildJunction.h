#ifndef JIAN_NUC3D_BUILDJUNCTION
#define JIAN_NUC3D_BUILDJUNCTION

#include "../geom/translate.h"
#include "HelixPar.h"
#include "../etl/ls.h"
#include "ScaffoldToAllAtom.h"

namespace jian {
namespace nuc3d {

template<typename ModelType = Model>
class BuildJunction {
public:
    using ResidueType = typename ModelType::ResidueType;

    std::mt19937 _rand_engine{11};
    std::uniform_real_distribution<double> _unif_distr{0, 1};

    nuc2d::loop *_loop;
    std::deque<ResidueType> _residues;

    ScaffoldToAllAtom<ModelType> scaffold_to_all_atom;

    BuildJunction() {}

    void init(nuc2d::loop *l) {
        _loop = l;
        _residues = init_residues();
    }

    void init(nuc2d::loop *l, const ModelType &model) {
        _loop = l;
        _residues = model.residues();
    }

    auto init_residues() {
        int num_branches = _loop->num_branches();
        int len = 4 * num_branches;
        auto dist_bound = mat::make_mat(len, len);
        DG::DihBoundType dih_bound;
        std::deque<decltype(_loop->head)> deq; int flag = 0;
        _loop->each([&](const auto &res, int index){
            deq.push_back(res);
            if (res->type == ')') flag++;
            if (flag == 2) {
                this->set_bound_by_helix(dist_bound, dih_bound, 
                    etl::ref(deq, -4)->num - 1, etl::ref(deq, -3)->num - 1, 
                    etl::ref(deq, -2)->num - 1, etl::ref(deq, -1)->num - 1);    
                for (int i = 0; i < 4; i++) deq.pop_back();
            }
        });
        return scaffold_to_residues(DG(dist_bound, dih_bound)());
    }

    template<typename MatType>
    std::deque<ResidueType> scaffold_to_residues(MatType &&mat) {}

    template<typename DistBound, typename DihBound>
    void set_bound_by_helix(DistBound &&dist_bound, DihBound &&dih_bound, int a1, int a2, int b1, int b2) {
        int len = mat::rows(dist_bound);
        int len_helix = a2 - a1 + 1;
        auto s = fpl::list<std::vector>(len_helix * 2, [&](int n){
            if (n < len_helix) return a1 + n; else return b1 + n - len_helix;});
        for (int n = 1; n < len_helix; n++) for (int i = 0; i < len_helix - 1 - n; i++) {
            int j = 2 * len_helix - 1 - i;
            mat::ref(dist_bound, s[i], s[i + n]) = mat::ref(dist_bound, s[i + n], s[i]) = HelixPar::dist_a(n);
            mat::ref(dist_bound, s[j - n], s[j]) = mat::ref(dist_bound, s[j], s[j - n]) = HelixPar::dist_a(n);
            mat::ref(dist_bound, s[i], s[j - n]) = mat::ref(dist_bound, s[j - n], s[i]) = HelixPar::dist_c(n);
            mat::ref(dist_bound, s[i + n], s[j]) = mat::ref(dist_bound, s[j], s[i + n]) = HelixPar::dist_d(n);
        }
        for (int i = 0; i < len_helix - 3; i++) {
            dih_bound(s[i], s[i + 1], s[i + 2], s[i + 3]) = HelixPar::dih_backbone;
            int j = i + len_helix;
            dih_bound(s[j], s[j + 1], s[j + 2], s[j + 3]) = HelixPar::dih_backbone;
        }
    }

    ModelType operator ()() {
        return build_junction();
    }

    ModelType build_junction() {
        auto ls = split();
        int num_branches = ls.size();
        int index = _unif_distr(_rand_engine) * num_branches;
        move(ls[index]);
        return pdb::residues_to_model(_residues);
    }

    auto split() {
        std::deque<std::deque<std::reference_wrapper<ResidueType>>> ls;
        if (_loop->is_open()) {
            _loop->each([](auto r, int index){});
        } else {
            _loop->each([](auto r, int index){});
        }
        return ls;
    }

    template<typename LS>
    void move(LS &&ls) {
        double flag = _unif_distr(_rand_engine);
        if (flag < 0.5) translate(ls);
        else rotate(ls);
    }

    template<typename LS>
    void translate(LS &&ls) {
        std::array<double, 3> s;
        for (int i = 0; i < 3; i++) s[i] = (_unif_distr(_rand_engine) - 0.5) * 2;
        for (auto && res : ls) geom::translate(res, s);
    }

    template<typename LS>
    void rotate(LS &&ls) {}
};

} // namespace nuc3d
} // namespace jian

#endif


