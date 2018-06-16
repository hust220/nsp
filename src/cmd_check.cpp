#include "nsp.hpp"
#include "pdb.hpp"

namespace jian {

REGISTER_NSP_COMPONENT(check) {
    RNA::check(mol_read_to<Model>(par.get("s")));
}

}

