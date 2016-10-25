#include <utility>
#include <cmath>
#include "THeadHairpin.hpp"

namespace jian {
namespace nuc3d {
namespace triple {

REGISTER_TRIPLE_MODULE_FACTORY("head_hairpin", THeadHairpin);

THeadHairpin::THeadHairpin(const Tuple &tuple, const Tuple &size) {
    int len = size[1];
    Tuple t = tuple;
    std::sort(t.begin(), t.end(), std::less<int>{});
    Frag frag;
        for (int i = 0; i <= t[0]; i++) {
            frag.push_back(i);
        }
        d_frags.push_back(std::move(frag));
        for (int i = t[1]; i <= t[2]; i++) {
            frag.push_back(i);
        }
        d_frags.push_back(std::move(frag));
        d_max_len = std::max(t[0], t[2] - t[1] - 1);
    d_indices = std::make_unique<Mati>(d_max_len, 3);
    for (int i = 0; i < d_max_len; i++) for (int j = 0; j < 3; j++) (*d_indices)(i, j) = -1;
    set_indices(0, d_frags[0].front(), d_frags[0].back() - 1);
    set_indices(1, d_frags[1].front() + 1, d_frags[1].back() - 1);
}

std::string THeadHairpin::type() const {
    return "head_hairpin";
}

} // namespace triple
}
} // namespace jian


