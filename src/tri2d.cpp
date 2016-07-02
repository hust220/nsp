#include "nsp.hpp"
#include <jian/nuc2d/PredTriSS.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(tri2d) {
    nuc2d::PredTriSS tri2d;
    int k = 4;
    par.set(k, "k");
    tri2d.run(par["seq"][0], k);
}

} // namespace jian

