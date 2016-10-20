#include "nsp.hpp"
#include <jian/nuc2d/ss_pred.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(ssp_fe) {
    std::string seq = par["seq"][0];
    std::string ss = ss_pred(seq);
    std::cout << seq << std::endl;
    std::cout << ss << std::endl;
}

} // namespace jian

