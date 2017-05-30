#include "nsp.hpp"
#include "pdb.hpp"
#include "tsdna.hpp"

BEGIN_JN

static void handle_tsdna_ssp(Str seq) {
    auto ls = tsdna_ssp(seq);
    for (auto && i : ls) {
        JN_OUT << i.ss << ' ' << i.score << std::endl;
    }
}

REGISTER_NSP_COMPONENT(tsdna) {
    auto g = par.getv("global");
    Str t = g[1];
    if (t == "ssp") {
        Str seq = g[2];
        handle_tsdna_ssp(seq);
    }
    else if (t == "tsp") {
        tsdna::THMC m;
        m.init(par);
        m.build_scaffold();
        m.run();
    }
}

END_JN

