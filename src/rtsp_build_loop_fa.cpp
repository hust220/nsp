#include "rtsp_build_loop_fa.hpp"

BEGIN_JN

Chain init_loop(Str seq, Str ss) {
    Int n = size(seq);
    Chain chain = bires_chain(seq[0], seq[1], ss[0], ss[1]);
    //mol_write(chain, "aa-1.pdb");
    for (Int i = 1; i < n-1; i++) {
        Chain c = bires_chain(seq[i], seq[i+1], ss[i], ss[i+1]);
        //mol_write(c, to_str("aa-", i+1, ".pdb"));
        auto sp = sp_bp(c[0], chain.back());
        for (auto && atom : c[1]) sp.apply(atom);
        chain.push_back(c[1]);
    }
    return chain;
}

Chain build_loop(Str seq, Str ss) {
    BuildLoopFA build_loop_fa;
    build_loop_fa.pred = init_loop(seq, ss);
    build_loop_fa.ss = ss;
//    build_loop_fa._mc_init_tempr = 1;
    build_loop_fa.run();
    return build_loop_fa.pred;
}

END_JN

