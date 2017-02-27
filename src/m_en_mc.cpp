#include "nsp.hpp"
#include <jian/geom.hpp>
#include <nsp/pdb.hpp>
#include <nsp/scoring/ScoreAa.hpp>
#include <nsp/scoring/Score.hpp>
#include <jian/utils/file.hpp>
#include <nsp/nuc3d/Format.hpp>
#include <jian/utils/log.hpp>

BEGIN_JN

namespace {

    struct En {
        Num len = 0, ang = 0, dih = 0, crash = 0, pairing = 0, stacking = 0;
    };

    Num en_total(const En &en) {
        return en.len + en.ang + en.dih + en.crash + en.pairing + en.stacking;
    }

    REGISTER_NSP_COMPONENT(en_mc) {
        Str filename;
        Chain chain;
        Int l;
        Num d;

        par.set(filename, "s");
        chain_read_model(chain, filename);
        l = size(chain);

        auto scorer = Score::fac_t::make("6p");
        scorer->init();

        En en;
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

        JN_OUT
            << "total: " << en_total(en)
            << ", len: " << en.len
            << ", ang: " << en.ang
            << ", dih: " << en.dih
            << ", crash: " << en.crash
            << ", pairing: " << en.pairing
            << ", stacking: " << en.stacking
            << std::endl;
    }
}

END_JN
















