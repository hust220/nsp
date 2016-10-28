#include "nsp.hpp"
#include <jian/lrsp/LRTSP.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(lrtsp) {
    lrsp::LRTSP tsp(par);
    tsp.pred();
}

} // namespace jian

