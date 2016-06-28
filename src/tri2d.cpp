#include "nsp.hpp"
#include <jian/nuc2d/PredTriSS.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(tri2d) {
    nuc2d::PredTriSS tri2d;
    tri2d.run(par["seq"][0]);
}

} // namespace jian

