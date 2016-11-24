#include "nsp.hpp"
#include <jian/nuc3d/BuildLoopRaw.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(build_loop) {
	int n = 1;
	std::string name = "aa";
	std::string seq = "AACCCCUU";
	std::string ss = "((....))";

	par.set(n, "n", "num");
	par.set(name, "name");
	par.set(seq, "seq");
	par.set(ss, "ss");

    BuildLoopRaw build_loop;
	build_loop.init(seq, ss);
    for (int i = 0; i < n; i++) {
        mol_write(build_loop(), name+'-'+JN_STR(i+1)+".pdb");
    }
}

} // namespace jian

