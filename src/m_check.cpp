#include "nsp.hpp"
#include <jian/pdb.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(check) {
    RNA::check(mol_read_to<Model>(par.get("s")));
}

END_JN

