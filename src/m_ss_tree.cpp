#include "nsp.hpp"
#include <jian/nuc2d/SSTree.hpp>
#include <jian/nuc2d/SSE.hpp>
#include <jian/utils/exception.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(ss_tree) {
    try {
		Str seq = "AAAACCCCUUUU";
		Str ss = "((((....))))";
        SSTree ss_tree;

		par.set(seq, "seq");
		par.set(ss, "ss");

        if (par.has("broken")) {
            ss_tree.make_b(seq, ss);
        } else {
            ss_tree.make(seq, ss);
        }
        ss_tree.head()->print_tree();
    } catch(const Error &e) {
        std::cout << e.what() << std::endl;
    }
}

END_JN

