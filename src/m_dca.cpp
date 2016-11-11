#include <jian/dca/Dca.hpp>
#include "nsp.hpp"

namespace jian {
namespace dca {

REGISTER_NSP_COMPONENT(dca) {
    int n = 1;
    float step = 1;
	std::string mol_type = "RNA";

    std::string out_file = par.get("o", "out");
    std::string fa_file = par.get("i", "in");
    std::string method = par.get("m", "method");

    par.set(n, "n");
    par.set(step, "step");
	par.set(mol_type, "t", "type");

    Dca *dca = FacDca::create(method);
    if (method == "mp") dca->set_step(step);
    dca->run(mol_type, fa_file, out_file, n-1);
    delete dca;
}

}
} // namespace jian

