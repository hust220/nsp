#include <utility>
#include <cmath>
#include "HeadHairpin.hpp"
#include "../../utils/Factory.hpp"

namespace jian {
namespace quadruple {

REGISTER_QUADRUPLE_MODULE_FACTORY("head_hairpin", HeadHairpin);

HeadHairpin::HeadHairpin(const Tuple &tuple, const Tuple &size) {
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
    for (int i = t[3]; i < len; i++) {
        frag.push_back(i);
    }
    d_frags.push_back(std::move(frag));
    d_max_len = std::max({t[0] + 1, t[2] - t[1] + 1, len - t[3]});
    d_indices = std::make_unique<Mat>(d_max_len, 4);
    for (int i = 0; i < d_max_len; i++) for (int j = 0; j < 4; j++) (*d_indices)(i, j) = -1;
    set_indices(0, d_frags[0].front(), d_frags[0].back());
    set_indices(1, d_frags[1].front(), d_frags[1].back());
    set_indices(2, d_frags[2].front(), d_frags[2].back());
}

} // namespace quadruple
} // namespace jian


