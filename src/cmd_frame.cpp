#include "nsp.hpp"
#include "pdb.hpp"
#include "geom.hpp"

BEGIN_JN

namespace {

REGISTER_NSP_COMPONENT(frame) {
    auto g = par.getv("global");
    auto && m = mol_read_to<Molecule>(g[1]);
    Int n = JN_INT(g[2])-1;

    JN_OUT << m[n] << Endl;
}

}

END_JN

