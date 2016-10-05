#include <jian/dca/Dca.hpp>
#include "nsp.hpp"

namespace jian {
namespace dca {

REGISTER_NSP_COMPONENT(dca) {
    int n = 1;
    par.set(n, "n");
    std::string out_file = par["out"][0];
    std::string fa_file = par["in"][0];
    std::string method = par["method"][0];

    Dca *dca = FacDca::create(method);
    dca->run(fa_file, out_file, n-1);
    delete dca;
}

}
} // namespace jian
