#include "nsp.hpp"
#include "rtsp_split.hpp"

BEGIN_JN

REGISTER_NSP_COMPONENT(split) {
    nuc3d::split(par);
}

END_JN

