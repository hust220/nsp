#include "nsp.hpp"
#include <jian/pdb.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(rna) {
    std::cout << mol_read_to<Model>(par.get("s", "pdb"), "RNA") << std::endl;
}

} // namespace jian

