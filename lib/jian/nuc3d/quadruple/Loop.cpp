#include "Loop.hpp"
#include "../../utils/Factory.hpp"

namespace jian {
namespace nuc3d {
namespace quadruple {

REGISTER_QUADRUPLE_MODULE_FACTORY("loop", Loop);

Loop::Loop(const Tuple &t1, const Tuple &t2) {
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
    d_max_len = std::accumulate(d_frags.begin(), d_frags.end(), 0, [](int n, auto &&frag){
        return std::max(n, int(frag.size()) - 2);
    });
    d_indices = std::make_unique<Mat>(d_max_len, 4);
    for (int i = 0; i < d_max_len; i++) for (int j = 0; j < 4; j++) (*d_indices)(i, j) = -1;
    for (int i = 0; i < 4; i++) {
        set_indices(i, p[i].first + 1, p[i].second - 1);
    }
}

std::string Loop::type() const {
    return "loop";
}

} // namespace quadruple
}
} // namespace jian


