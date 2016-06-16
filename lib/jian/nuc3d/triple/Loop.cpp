#include "Loop.hpp"
#include "../../utils/Factory.hpp"

namespace jian {
namespace triple {

REGISTER_TRIPLE_MODULE_FACTORY("loop", Loop);

Loop::Loop(const Tuple &t1, const Tuple &t2) {
    Tuple s1 = t1, s2 = t2;
    std::sort(s1.begin(), s1.end(), std::less<int>{});
    std::sort(s2.begin(), s2.end(), std::less<int>{});
    Frag frag;
    std::pair<int, int> p[3];
    for (int i = 0; i < 3; i++) {
        p[i] = std::minmax(s1[i], s2[i]);
        for (int j = p[i].first; j <= p[i].second; j++) {
            frag.push_back(j);
        }
        d_frags.push_back(std::move(frag));
    }
    d_max_len = std::accumulate(d_frags.begin(), d_frags.end(), 0, [](int n, auto &&frag){
        return std::max(n, int(frag.size()));
    });
    d_indices = std::make_unique<Mat>(d_max_len, 3);
    for (int i = 0; i < d_max_len; i++) for (int j = 0; j < 3; j++) (*d_indices)(i, j) = -1;
    for (int i = 0; i < 3; i++) {
        set_indices(i, p[i].first, p[i].second);
    }
}

std::string Loop::type() const {
    return "loop";
}

} // namespace triple
} // namespace jian


