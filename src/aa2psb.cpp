#include <sstream>
#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/cg.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(aa2psb) {
    std::cout << CGpsb::chain(read_model_to_chain(par.get("s"))) << std::endl;
}

} // namespace jian
















