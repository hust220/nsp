#include <string>
#include "nsp.hpp"
#include <jnbio/nuc3d/transform.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(transform) {
    auto && m = mol_read_to<Model>(par.get("s"));
    std::cout << transform(m, par["seq"][0], par["type"][0]) << std::endl;
}

END_JN

