#include "nsp.hpp"
#include "rss_ss_pred.hpp"

namespace jian {

REGISTER_NSP_COMPONENT(ssp_fe) {
    S seq = par.get("seq");
    S ss = ss_pred(seq);
    JN_OUT << seq << Endl;
    JN_OUT << ss  << Endl;
}

}

