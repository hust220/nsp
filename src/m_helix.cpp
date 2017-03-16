#include "nsp.hpp"
#include <nsp/rtsp/build_helix.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(helix) {
    auto globals = par.getv("global");
    Str m = globals[1];
    if (m == "build") {
        Str seq = globals[2];
        JN_OUT << build_helix(seq);
    }
}

END_JN

