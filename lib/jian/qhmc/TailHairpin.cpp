#include <utility>
#include <cmath>
#include "TailHairpin.hpp"

namespace jian {
namespace qhmc {

REGISTER_QUADRUPLE_MODULE_FACTORY("tail_hairpin", TailHairpin);

TailHairpin::TailHairpin(const Tuple &tuple, const Tuple &size) {
    int len = size[1];
    Tuple t = tuple;
    std::sort(t.begin(), t.end(), std::less<int>{});
    Frag frag;
    for (int i = t[0]; i <= t[1]; i++) {
        frag.push_back(i);
    }
    d_frags.push_back(std::move(frag));
    for (int i = t[2]; i <= t[3]; i++) {
        frag.push_back(i);
    }
    d_frags.push_back(std::move(frag));
    d_max_len = std::max(t[1] - t[0] - 1, t[3] - t[2] - 1);
    d_indices = std::make_unique<Mat>(d_max_len, 4);
    for (int i = 0; i < d_max_len; i++) for (int j = 0; j < 4; j++) (*d_indices)(i, j) = -1;
    set_indices(0, d_frags[0].front() + 1, d_frags[0].back() - 1);
    set_indices(1, d_frags[1].front() + 1, d_frags[1].back() - 1);
}

std::string TailHairpin::type() const {
    return "tail_hairpin";
}

} // namespace quadruple
}


