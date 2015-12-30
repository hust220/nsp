#ifndef JIAN_DG_DIHBOUNDS
#define JIAN_DG_DIHBOUNDS

#include "../util/std.h"
#include "../util/MLib.h"

namespace jian {
namespace dg {

template<typename LS>
class hash {
public:
    unsigned int operator ()(const LS &vec) const {
        return fold([](unsigned int result, int i){return result ^ (std::hash<int>{}(i) << 1);}, 0, vec);
    }
};

template<typename LS>
class equal_to {
public:
    bool operator ()(const LS &vec1, const LS &vec2) const {
        return vec1.size() == vec2.size() and not exists([](auto i, auto j){return i != j;}, vec1, vec2);
    }
};

template<typename LS = std::vector<int>>
class DihBound : public std::unordered_map<LS, std::pair<double, double>, hash<LS>, equal_to<LS>> {
public:

    std::pair<double, double> &operator ()(int a, int b, int c, int d) {
        return (*this)[{a, b, c, d}];
    }

    bool exists(int a, int b, int c, int d) const {
        return (*this).count({a, b, c, d});
    }
};

template<typename LS>
std::ostream &operator <<(std::ostream &out, const DihBound<LS> &dih_bound) {
    for (auto && pair : dih_bound) {
        for (auto i : pair.first) std::cout << i << ' ';
        std::cout << ':' << pair.second.first << ' ' << pair.second.second << std::endl;
    }
}

} // namespace dg
} // namespace jian

#endif

