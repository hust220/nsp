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

END_JN

