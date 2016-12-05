#include "nsp.hpp"
#include <jian/lrsp/LRTSP.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(lrtsp) {
    lrsp::LRTSP tsp(par);
    tsp.pred();
}

END_JN

