#include "nsp.hpp"
#include <jian/nuc2d/SSTree.hpp>
#include <jian/nuc2d/loop.hpp>
#include <jian/utils/exception.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(ss_tree) {
    try {
		str_t seq = "AAAACCCCUUUU";
		str_t ss = "((((....))))";
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

} // namespace jian

