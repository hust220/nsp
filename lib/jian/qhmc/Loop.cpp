#include "Loop.hpp"

namespace jian {
namespace qhmc {

REGISTER_QUADRUPLE_MODULE_FACTORY("loop", Loop);

Loop::Loop(const Tuple &t1, const Tuple &t2, int len) {
    Tuple s1 = t1, s2 = t2;
    Frag frag;
    std::pair<int, int> p[4];
    for (int i = 0; i < 4; i++) {
        p[i] = std::minmax(s1[i], s2[i]);
        for (int j = p[i].first+1; j < p[i].second; j++) {
            frag.push_back(j);
        }
        d_frags.push_back(std::move(frag));
    }
    d_max_len = std::accumulate(d_frags.begin(), d_frags.end(), 0, [](int n, auto &&frag){
        return std::max(n, int(frag.size()));
    });
    d_indices = Mat::Constant(d_max_len, 4, -1);
    for (int i = 0; i < 4; i++) {
        set_indices(i, d_frags[i], t1[i] <= t2[i]);
    }
}

std::string Loop::type() const {
    return "loop";
}

} // namespace quadruple
}


