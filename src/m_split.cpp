#include "nsp.hpp"
#include <nsp/nuc3d/split.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(split) {
    nuc3d::split(par);
}

END_JN
