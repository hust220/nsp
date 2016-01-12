#ifndef JIAN_NUC3D_TRIASS_BUILDTRILOOP
#define JIAN_NUC3D_TRIASS_BUILDTRILOOP

#include "../../dg/DG.h"

namespace jian {
namespace nuc3d {
namespace triass {

template<typename ModelType>
class BuildTriLoop {
public:
    template<typename LS>
    ModelType operator ()(const LS &ls, int primary_index, int secondary_index) {
        return build_tri_loop(ls, primary_index, secondary_index);
    }

    template<typename LS>
    ModelType build_tri_loop(const LS &ls, int primary_index, int secondary_index) {
        DG::DistBoundType dist_bound; DG::DihBoundType dih_bound;
        std::tie(dist_bound, dih_bound) = bound_constraints(ls, primary_index, secondary_index);
        DG dg(dist_bound, dih_bound);
        auto scaffold = dg();
        std::cout << scaffold << std::endl;
        return all_atom(scaffold);
    }

    template<typename LS>
    std::pair<DG::DistBoundType, DG::DihBoundType> bound_constraints(const LS &ls, int primary_index, int secondary_index) {
        int len = fold([](double sum, int n){return sum + n;}, 0, ls) + 6;
        auto dist_bound = mat::make_mat<DG::DistBoundType>(len, len);
        DG::DihBoundType dih_bound;

        init_bound(dist_bound, len);
        set_chain(dist_bound, dih_bound, let([&ls]{
            int index = 0; return map([&index](int n){index += n + 2; return range(index - n - 2, n + 2);}, ls);
        }));
        set_triple(dist_bound, dih_bound, let([&ls, &len]{
            auto l = map<std::vector>([](int i){return i;}, ls); 
            return std::vector<std::vector<int>>{{0, l[0] + 2, l[0] + l[1] + 4}, {l[0] + 1, l[0] + l[1] + 3, len - 1}};
        }), primary_index, secondary_index);

        return {dist_bound, dih_bound};
    }

    void init_bound(DG::DistBoundType &dist_bound, int len) {
        for (int i = 0; i < len; i++) {
            for (int j = i + 1; j < len; j++) {
                mat::ref(dist_bound, i, j) = 99;
                mat::ref(dist_bound, j, i) = 6.1;
            }
        }
    }

    template<typename LS>
    void set_chain(DG::DistBoundType &dist_bound, DG::DihBoundType &dih_bound, const LS &ls) {
        for (auto && v : ls) for (int i = 0; i < v.size() - 1; i++) {
            dist_bound(v[i], v[i + 1]) = dist_bound(v[i + 1], v[i]) = 6.1;
            if (i < v.size() - 2) dist_bound(v[i], v[i + 2]) = dist_bound(v[i + 2], v[i]) = 11;
        }
    }

    template<typename LS>
    void set_triple(DG::DistBoundType &dist_bound, DG::DihBoundType &dih_bound, const LS &ls, int primary_index, int secondary_index) {
        int i1 = primary_index, i2 = secondary_index, i3 = 3 - i1 - i2;
        auto at = [&ls](int i, int j){
            return (j % 2 == i ? ls[i][j] + 1 : ls[i][j] - 1);
        };
        dist_bound(ls[0][i1], ls[0][i2]) = dist_bound(ls[0][i2], ls[0][i1]) = 15.1;
        dist_bound(ls[0][i1], at(0, i2)) = dist_bound(at(0, i2), ls[0][i1]) = 18.0;
        dist_bound(ls[0][i2], at(0, i1)) = dist_bound(at(0, i1), ls[0][i2]) = 12.4;
        dih_bound(ls[0][i1], ls[0][i3], at(0, i3), ls[0][i2]) = {-0.31562, -0.31562};

        dist_bound(ls[0][i1], ls[0][i3]) = dist_bound(ls[0][i3], ls[0][i1]) = 11.1;
        dist_bound(ls[0][i1], at(0, i3)) = dist_bound(at(0, i3), ls[0][i1]) = 7.3;
        dist_bound(ls[0][i3], at(0, i1)) = dist_bound(at(0, i1), ls[0][i3]) = 14.8;
        dih_bound(ls[0][i1], ls[0][i2], at(0, i2), ls[0][i3]) = {1.35937, 1.35937};

        dist_bound(ls[0][i2], ls[0][i3]) = dist_bound(ls[0][i3], ls[0][i2]) = 18.1;
        dist_bound(ls[0][i2], at(0, i3)) = dist_bound(at(0, i3), ls[0][i2]) = 18.4;
        dist_bound(ls[0][i3], at(0, i2)) = dist_bound(at(0, i2), ls[0][i3]) = 18.2;
        dih_bound(ls[0][i2], ls[0][i1], at(0, i1), ls[0][i3]) = {-0.869972, -0.869972};

        dist_bound(ls[1][i1], ls[1][i2]) = dist_bound(ls[1][i2], ls[1][i1]) = 15.1;
        dist_bound(ls[1][i1], at(1, i2)) = dist_bound(at(1, i2), ls[1][i1]) = 12.4;
        dist_bound(ls[1][i2], at(1, i1)) = dist_bound(at(1, i1), ls[1][i2]) = 18.0;
        dih_bound(ls[1][i1], ls[1][i3], at(1, i3), ls[1][i2]) = {0.900876, 0.900876};

        dist_bound(ls[1][i1], ls[1][i3]) = dist_bound(ls[1][i3], ls[1][i1]) = 11.1;
        dist_bound(ls[1][i1], at(1, i3)) = dist_bound(at(1, i3), ls[1][i1]) = 14.8;
        dist_bound(ls[1][i3], at(1, i1)) = dist_bound(at(1, i1), ls[1][i3]) = 7.3;
        dih_bound(ls[1][i1], ls[1][i2], at(1, i2), ls[1][i3]) = {-1.37275, -1.37275};

        dist_bound(ls[1][i2], ls[1][i3]) = dist_bound(ls[1][i3], ls[1][i2]) = 18.1;
        dist_bound(ls[1][i2], at(1, i3)) = dist_bound(at(1, i3), ls[1][i2]) = 18.2;
        dist_bound(ls[1][i3], at(1, i2)) = dist_bound(at(1, i2), ls[1][i3]) = 18.4;
        dih_bound(ls[1][i2], ls[1][i1], at(1, i1), ls[1][i3]) = {0.296515, 0.296515};
    }

    ModelType all_atom(const MatrixXf &dist) {
        return ModelType();
    }
};


} // namespace triass
} // namespace nuc3d
} // namespace jian

#endif

