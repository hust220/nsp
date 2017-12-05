#include <utility>
#include <cmath>
#include <numeric>
#include "g4dna_head_hairpin.hpp"

BEGIN_JN
namespace qhmc {

REGISTER_QUADRUPLE_MODULE_FACTORY("head_hairpin", HeadHairpin);

HeadHairpin::HeadHairpin(const Tuple &head, const Tuple &tail, int len) {
    int i;
    Frag frag;

    for (i = 0; i < head[0]; i++) {
        frag.push_back(i);
    }
    d_frags.push_back(std::move(frag));

    if (head[1] >= tail[1]) {
        if (head[2] < tail[2]) {
            for (i = head[1]+1; i < head[2]; i++) {
                frag.push_back(i);
            }
        } else {
            for (i = head[1]+1; i < tail[2]; i++) {
                frag.push_back(i);
            }
        }
    }
    d_frags.push_back(std::move(frag));

    if (head[2] >= tail[2]) {
        if (head[3] < tail[3]) {
            for (i = head[2]+1; i < head[3]; i++) {
                frag.push_back(i);
            }
        } else {
            for (i = head[2]+1; i < tail[3]; i++) {
                frag.push_back(i);
            }
        }
    }
    d_frags.push_back(std::move(frag));

    if (head[3] >= tail[3]) {
        for (i = head[3]+1; i < len; i++) {
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

S HeadHairpin::type() const {
    return "head_hairpin";
}

}
}


