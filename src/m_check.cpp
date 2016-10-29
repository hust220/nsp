#include "nsp.hpp"
#include <jian/pdb.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(check) {
    RNA::check(mol_read_to<Model>(par.get("s")));
}

} // namespace jian

