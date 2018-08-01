#include "nsp.hpp"
#include "rtsp_build_loop_dg.hpp"

namespace jian {

REGISTER_NSP_COMPONENT(build_loop_dg) {
    Str seq = par.get("seq");
    Str ss = par.get("ss");
    JN_OUT << *build_chain_dg(seq, ss) << std::endl;
}

}

