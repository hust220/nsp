#include "nsp.hpp"
#include <jian/pdb.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(check) {
    Model m(par["pdb"][0]);
    RNA::check(m);
}

} // namespace jian

