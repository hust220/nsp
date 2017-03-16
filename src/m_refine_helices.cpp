#include "nsp.hpp"
#include <nsp/pdb.hpp>

BEGIN_JN

namespace {

    void refine_helices(Molecule &m, const Str &ss) {
    }

    REGISTER_NSP_COMPONENT(refine_helices) {
        Str ss = par.get("ss");
        Str filename = par.get("s");

        Molecule m;
        mol_read(m, filename);

        refine_helices(m, ss);

        JN_OUT << m << std::endl;
    }

}

END_JN

