#include <string>
#include "nsp.hpp"
#include <jian/nuc3d/transform.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(transform) {
    Model m(par["pdb"][0]);
    std::cout << transform(m, par["seq"][0], par["type"][0]) << std::endl;
}

} // namespace jian

