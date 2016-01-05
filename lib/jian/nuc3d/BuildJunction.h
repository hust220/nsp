#ifndef JIAN_NUC3D_BUILDJUNCTION
#define JIAN_NUC3D_BUILDJUNCTION

#include "../geom/translate.h"

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

    BuildJunction() {}

    void init(nuc2d::loop *l) {
        _loop = l;
        _residues = init_residues();
    }

    void init(nuc2d::loop *l, const ModelType &model) {
        _loop = l;
        _residues = model.residues();
    }

    std::deque<ResidueType> init_residues() {
        return std::deque<ResidueType>();
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


