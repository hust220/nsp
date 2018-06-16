#pragma once

#include "score_par_bp.hpp"
#include "dca.hpp"

namespace jian {

template<typename _Chain>
Str get_ss(const _Chain &chain) {
    ParBp pb;
    int i, j, l;
    ::jian::dca::pairs_t pairs;

    l = size(chain);
    for (i = 0; i < l; i++) {
        for (j = i + 1; j < l; j++) {
            if (is_bp(chain[i], chain[j])) {
                pairs.push_back({ i, j });
            }
        }
    }
    //::jian::dca::print_pairs(pairs);
    return ::jian::dca::pairs_to_ss(pairs, chain.size());
}

}

