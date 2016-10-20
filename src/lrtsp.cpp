#include "nsp.hpp"
#include <jian/lrsp/TSP.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(lrtsp) {
    lrsp::TSP tsp(par);
    tsp.pred();
}

} // namespace jian

