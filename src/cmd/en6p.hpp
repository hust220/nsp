#pragma once

#include "score.hpp"

namespace jian {

struct En6p {
    Num len = 0, ang = 0, dih = 0, crash = 0, pairing = 0, stacking = 0, rg = 0;
};

inline Num total(const En6p &en) {
    return en.len + 5*en.ang + 5*en.dih + en.crash + en.pairing + en.stacking + 0.02*en.rg;
}

template<typename _Chain>
En6p en6p_chain(const _Chain &c) {
    auto scorer = Score::fac_t::make("6p");
    scorer->init();

    auto cg = CG::fac_t::make("6p");

    Chain chain;
    for (auto && r : c) chain.push_back(cg->to_cg(r));

    Int l = size(chain);
    Num d;

    En6p en;
    en.rg = en_rg_6p(chain);
    for (Int i = 0; i < l; i++) {
        en.len += scorer->en_len(chain, i);
        en.ang += scorer->en_ang(chain, i);
        en.dih += scorer->en_dih(chain, i);
        for (Int j = i + 1; j < l; j++) {
            if (geom::distance(chain[i][2], chain[j][2]) < 20) {
                en.crash += scorer->en_crash(chain[i], chain[j]);
                scorer->en_bp(chain[i], chain[j]);
                en.pairing += scorer->m_en_pairing;
                en.stacking += scorer->m_en_stacking;
            }
        }
    }

    return en;

}


}

