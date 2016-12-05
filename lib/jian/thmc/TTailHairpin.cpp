#include <utility>
#include <cmath>
#include "TTailHairpin.hpp"

BEGIN_JN
namespace nuc3d {
namespace triple {

REGISTER_TRIPLE_MODULE_FACTORY("tail_hairpin", TTailHairpin);

TTailHairpin::TTailHairpin(const Tuple &tuple, const Tuple &size) {
    int len = size[1];
    Tuple t = tuple;
    std::sort(t.begin(), t.end(), std::less<int>{});
    Frag frag;
        for (int i = t[0]; i <= t[1]; i++) {
            frag.push_back(i);
        }
        d_frags.push_back(std::move(frag));
        for (int i = t[2]; i < len; i++) {
            frag.push_back(i);
        }
        d_frags.push_back(std::move(frag));
        d_max_len = std::max(t[1] - t[0] - 1, len - t[2] - 1);
    d_indices = std::make_unique<Mati>(d_max_len, 3);
    for (int i = 0; i < d_max_len; i++) for (int j = 0; j < 3; j++) (*d_indices)(i, j) = -1;
    set_indices(0, d_frags[0].front() + 1, d_frags[0].back() - 1);
    set_indices(1, d_frags[1].front() + 1, len - 1);
}

S TTailHairpin::type() const {
    return "tail_hairpin";
}

} // namespace triple
}
END_JN


