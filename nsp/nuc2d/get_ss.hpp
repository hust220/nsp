#pragma once

#include "../scoring/ParBp.hpp"
#include "../dca.hpp"

BEGIN_JN

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

END_JN

