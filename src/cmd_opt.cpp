#include "nsp.hpp"
#include "rna_mc.hpp"
#include "tsdna.hpp"
#include "g4dna.hpp"

BEGIN_JN

REGISTER_NSP_COMPONENT(opt) {
    auto g = par.getv("global");
    auto dhmc = DHMC::make(par);
    if (size(g) > 1 && g[1] == "score") {
        dhmc->score();
    }
    else {
        dhmc->predict();
    }
}

ALIAS_NSP_COMPONENT(opt, dhmc);

REGISTER_NSP_COMPONENT(thmc) {
    tsdna::THMC tri;
    tri.init(par);
    tri.predict();
}

REGISTER_NSP_COMPONENT(qhmc) {
    qhmc::QHMC qua;
    qua.init(par);
    qua.predict();
}

END_JN

