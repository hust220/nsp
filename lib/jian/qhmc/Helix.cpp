#include "Helix.hpp"

namespace jian {
namespace nuc3d {
namespace quadruple {

REGISTER_QUADRUPLE_MODULE_FACTORY("helix", Helix);

Helix::Helix(const Tuple &t1, const Tuple &t2) {
    Tuple s1 = t1, s2 = t2;
    Frag frag;
    std::pair<int, int> p[4];
    for (int i = 0; i < 4; i++) {
        p[i] = std::minmax(s1[i], s2[i]);
        for (int j = p[i].first; j <= p[i].second; j++) {
            frag.push_back(j);
        }
        d_frags.push_back(std::move(frag));
    }
    d_max_len = p[0].second - p[0].first + 1;
    d_indices = std::make_unique<Mat>(d_max_len, 4);
    for (int i = 0; i < d_max_len; i++) for (int j = 0; j < 4; j++) (*d_indices)(i, j) = -1;
    for (int i = 0; i < 4; i++) {
        set_indices(i, p[i].first, p[i].second);
    }
}

std::string Helix::type() const {
    return "helix";
}

} // namespace quadruple
}
} // namespace jian


