#include "nsp.hpp"
#include <nsp/nuc2d/ss_pred.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(ssp_fe) {
    S seq = par.get("seq");
    S ss = ss_pred(seq);
    JN_OUT << seq << Endl;
    JN_OUT << ss  << Endl;
}

END_JN

