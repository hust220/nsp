#include <utility>
#include <cmath>
#include "TailHairpin.hpp"
#include <numeric>

BEGIN_JN
namespace qhmc {

REGISTER_QUADRUPLE_MODULE_FACTORY("tail_hairpin", TailHairpin);

TailHairpin::TailHairpin(const Tuple &head, const Tuple &tail, int len) {
    int i;
    Frag frag;

    if (head[1] >= tail[1]) {
        for (i = tail[0]+1; i < tail[1]; i++) {
            frag.push_back(i);
        }
    } else {
        for (i = tail[0]+1; i < head[1]; i++) {
            frag.push_back(i);
        }
    }
    d_frags.push_back(std::move(frag));

    if (head[1] < tail[1]) {
        if (head[2] < tail[2]) {
            for (i = tail[1]+1; i < head[2]; i++) {
                frag.push_back(i);
            }
        } else {
            for (i = tail[1]+1; i < tail[2]; i++) {
                frag.push_back(i);
            }
        }
    }
    d_frags.push_back(std::move(frag));

    if (head[2] < tail[2]) {
        if (head[3] < tail[3]) {
            for (i = tail[2]+1; i < head[3]; i++) {
                frag.push_back(i);
            }
        } else {
            for (i = tail[2]+1; i < tail[3]; i++) {
                frag.push_back(i);
            }
        }
    }
    d_frags.push_back(std::move(frag));

    if (head[3] < tail[3]) {
        for (i = tail[3]; i < len; i++) {
            frag.push_back(i);
        }
    }
    d_frags.push_back(std::move(frag));

    d_max_len = std::accumulate(d_frags.begin(), d_frags.end(), 0, [&](int l, auto && frag){return std::max(l, int(frag.size()));});

    // Set indices
    d_indices = Mati::Constant(d_max_len, 4, -1);
    for (int i = 0; i < 4; i++) {
        set_indices(i, d_frags[i], head[i] <= tail[i]);
    }
}

S TailHairpin::type() const {
    return "tail_hairpin";
}

} // namespace quadruple
}


