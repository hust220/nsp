#include "nsp.hpp"
#include <jian/nuc3d/Assemble.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(assemble) {
    nuc3d::Assemble ass(par);
    ass.predict();
}

END_JN

