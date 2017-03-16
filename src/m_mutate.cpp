#include <string>
#include "nsp.hpp"
#include <nsp/rtsp/mutate.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(mutate) {
    auto global = par.getv("global");
    JN_OUT << mutate(mol_read_to<Chain>(global[1]), global[2], "RNA");
}

END_JN

