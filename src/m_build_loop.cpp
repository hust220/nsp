#include "nsp.hpp"
#include <nsp/nuc3d/BuildLoopRaw.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(build_loop) {
	int n = 1;
	Str name = "aa";
	Str seq = "AACCCCUU";
	Str ss = "((....))";

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

END_JN

