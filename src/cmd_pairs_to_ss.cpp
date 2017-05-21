#include "rssp.hpp"
#include "dca.hpp"
#include "log.hpp"
#include "nsp.hpp"
#include "jian.hpp"

BEGIN_JN

REGISTER_NSP_COMPONENT(pairs_to_ss) {
    Int l;
    par.set(l, "l");

    Str di_file = par.get("p");

    auto && pairs = dca::pairs_from_file(di_file, l);

    JN_OUT << dca::pairs_to_ss(pairs, l) << std::endl;

}

END_JN

