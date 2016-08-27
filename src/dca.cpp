#include <jian/dca.hpp>
#include "nsp.hpp"

namespace jian {

REGISTER_NSP_COMPONENT(dca) {
    int n = 0;
    par.set(n, "n");
    dca::result_t rt;
    dca::analyze(par.get("fasta"), n-1, rt);
    dca::print_dis(rt.dis);
}

} // namespace jian

