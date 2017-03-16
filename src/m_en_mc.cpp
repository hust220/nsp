#include "nsp.hpp"
#include <jian/geom.hpp>
#include <nsp/pdb.hpp>
#include <nsp/scoring/ScoreAa.hpp>
#include <nsp/scoring/Score.hpp>
#include <jian/utils/file.hpp>
#include <nsp/nuc3d/Format.hpp>
#include <jian/utils/log.hpp>
#include "En6p.hpp"

BEGIN_JN

namespace {

    REGISTER_NSP_COMPONENT(en_mc) {
        Str filename;
        par.set(filename, "s");

        Chain chain;
        chain_read_model(chain, filename);

        auto en = en6p_chain(chain);

        JN_OUT
            << "total: " << total(en)
            << ", len: " << en.len
            << ", ang: " << en.ang
            << ", dih: " << en.dih
            << ", crash: " << en.crash
            << ", pairing: " << en.pairing
            << ", stacking: " << en.stacking
            << ", rg: " << en.rg
            << std::endl;
    }
}

END_JN
















