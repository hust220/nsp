#include "nsp.hpp"
#include "rtsp.hpp"

BEGIN_JN

REGISTER_NSP_COMPONENT(frag) {
    auto && globals = par.getv("global");
    for (auto && i : globals) JN_OUT << i << ' '; JN_OUT << std::endl;
    if (globals[1] == "build") {
        JN_OUT << build_loop(globals[2], globals[3]);
    }
    else if (globals[1] == "sample") {
        Chain strand;
        mol_read(strand, globals[2]);
        sample_loop(strand, globals[3]);
        JN_OUT << strand;
    }
}

END_JN

