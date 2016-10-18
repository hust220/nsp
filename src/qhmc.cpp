#include "nsp.hpp"
#include <jian/qhmc/QHMC.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(qhmcpsb) {
    qhmc::QHMC<CGpsb> qua;
    qua.init(par);
    qua.run_job();
}

REGISTER_NSP_COMPONENT(qhmc1p) {
    qhmc::QHMC<CG1p> qua;
    qua.init(par);
    qua.run_job();
}

} // namespace jian

