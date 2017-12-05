#include "nsp.hpp"
#include "rtsp.hpp"

BEGIN_JN

REGISTER_NSP_COMPONENT(build_loop_fa) {
    BuildLoopFA build_loop_fa;
    Str seq = par.get("seq");
    Str ss = par.get("ss");

    build_loop_fa.pred = init_loop(seq, ss);
    build_loop_fa.ss = ss;
    build_loop_fa.show_log = true;
    build_loop_fa.traj_name = par.get("traj");
    build_loop_fa.run();

    JN_OUT << build_loop_fa.pred;

//    auto && g = par.getv("global");
//    if (g[1] == "build") {
//        Str seq = g[2];
//        Str ss = g[3];
//        JN_OUT << build_loop(seq, ss);
//    }
//    else if (g[1] == "sample") {
//        Str struct_name = g[2];
//        Str ss = g[3];
//
//        Chain strand;
//        mol_read(strand, struct_name);
//        sample_loop(strand, ss);
//        JN_OUT << strand;
//    }
}

END_JN

