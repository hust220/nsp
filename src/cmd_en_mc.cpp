#include "nsp.hpp"
#include "geom.hpp"
#include "pdb.hpp"
#include "score_aa.hpp"
#include "score.hpp"
#include "file.hpp"
#include "rtsp_format.hpp"
#include "log.hpp"
#include "cmd_en6p.hpp"

namespace jian {

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

}
















