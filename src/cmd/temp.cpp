#include "nsp.hpp"
#include "rna_mc.hpp"

namespace jian {

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

}

