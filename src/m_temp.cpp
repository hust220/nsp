#include "nsp.hpp"
#include <nsp/dhmc/DHMC.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(temp) {
    Str nat;
    Str prefix;
    Num rmsd;

    nat = par.get("nat", "n");
    rmsd = lexical_cast<Num>(par.get("rmsd", "r"));
    prefix = par.get("prefix", "p");

    auto dhmc = DHMC::make(par);
    dhmc->predict();
}

END_JN

