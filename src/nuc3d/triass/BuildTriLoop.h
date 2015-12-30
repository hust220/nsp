#ifndef JIAN_NUC3D_TRIASS_BUILDTRILOOP
#define JIAN_NUC3D_TRIASS_BUILDTRILOOP

namespace jian {
namespace nuc3d {
namespace triass {

template<typename ModelType>
class BuildTriLoop {
public:
    template<template<typename...> class ListType>
    ModelType operator ()(const ListType<int> &ls, int primary_index, int secondary_index) {
        return build_tri_loop(ls, primary_index, secondary_index);
    }

    template<template<typename...> class ListType>
    ModelType build_tri_loop(const ListType<int> &ls, int primary_index, int secondary_index) {
        DG::DistType dist_bound; DG::DihType dih_bound;
        std::tie(dist_bound, dih_bound) = bound_constraints(ls, primary_index, secondary_index);
        DG dg(dist_bound, dih_bound);
        auto scaffold = dg();
        return all_atom(scaffold);
    }

    template<template<typename...> class ListType>
    std::pair<DG::DistType, DG::DihType> bound_constraints(const ListType<int> &ls, int primary_index, int secondary_index) {
        int len = fold([](double sum, int n){return sum + n;}, 0, ls) + 6;
        auto dist_bound = mat::make_mat<DG::DistType>(len, len);
        set_bound_chain(dist_bound, let([&ls]{
            int index = 0; return map([&index](int n){index += n + 2; return range(index - n - 2, n + 2);}, ls);
        }));
        set_bound_triple(dist_bound, let([&ls, &len]{
            auto l = map<std::vector>([](int i){return i;}, ls); 
            return std::vector<std::vector<int>>{{0, l[0] + 2, l[0] + l[1] + 4}, {l[0] + 1, l[0] + l[1] + 3, len - 1}};
        }), primary_index, secondary_index);
        DG::DihType dih_bound;
        return {dist_bound, dih_bound};
    }

    template<typename LS>
    void set_bound_chain(DG::DistType &dist_bound, const LS &ls) {
        for (auto && v : ls) for (int i = 0; i < v.size() - 1; i++) {
            dist_bound(v[i], v[i + 1]) = dist_bound(v[i + 1], v[i]) = 6.1;
            if (i < v.size() - 2) dist_bound(v[i], v[i + 2]) = dist_bound(v[i + 2], v[i]) = 11;
        }
    }

    template<typename LS>
    void set_bound_triple(DG::DistType &dist_bound, const LS &ls, int primary_index, int secondary_index) {
        int i1 = primary_index, i2 = secondary_index, i3 = 3 - i1 - i2;
        for (auto i : {0, 1}) dist_bound(ls[i][i1], ls[i][i2]) = dist_bound(ls[i][i2], ls[i][i1]) = 15.1;
        for (auto i : {0, 1}) dist_bound(ls[i][i1], ls[i][i3]) = dist_bound(ls[i][i3], ls[i][i1]) = 15.1;
        for (auto i : {0, 1}) dist_bound(ls[i][i3], ls[i][i2]) = dist_bound(ls[i][i2], ls[i][i3]) = 18.1;
    }

    ModelType all_atom(const MatrixXf &dist) {
        return ModelType();
    }
};


} // namespace triass
} // namespace nuc3d
} // namespace jian

#endif

