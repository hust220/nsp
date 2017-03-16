#pragma once

#include <nsp/scoring/Score.hpp>

BEGIN_JN

struct En6p {
    Num len = 0, ang = 0, dih = 0, crash = 0, pairing = 0, stacking = 0, rg = 0;
};

inline Num total(const En6p &en) {
    return en.len + 5*en.ang + 5*en.dih + en.crash + en.pairing + en.stacking + 0.02*en.rg;
}

template<typename _Chain>
En6p en6p_chain(const _Chain &c) {
    Chain chain;
    for (auto && r : c) chain.push_back(r);

    Int l = size(chain);
    Num d;

    auto scorer = Score::fac_t::make("6p");
    scorer->init();

    En6p en;
    en.rg = en_rg(chain);
    for (Int i = 0; i < l; i++) {
        en.len += scorer->en_len(chain, i);
        en.ang += scorer->en_ang(chain, i);
        en.dih += scorer->en_dih(chain, i);
        for (Int j = i + 1; j < l; j++) {
            en.crash += scorer->en_crash(chain[i], chain[j]);
            scorer->en_bp(chain[i], chain[j]);
            en.pairing += scorer->m_en_pairing;
            en.stacking += scorer->m_en_stacking;
        }
    }

    return en;

}


END_JN

