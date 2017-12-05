#include "nsp.hpp"
#include "rss_ss_tree.hpp"
#include "rss_sse.hpp"
#include "exception.hpp"

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
        std::cout << ss_tree << std::endl;
        //ss_tree.head()->print_tree();
    } catch(const Error &e) {
        std::cout << e.what() << std::endl;
    }
}

END_JN

