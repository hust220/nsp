#include <sstream>
#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/cg.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(aa2psb) {
    std::cout << CGpsb::chain(residues_from_file(par["pdb"][0])) << std::endl;
}

} // namespace jian
















