#include "nsp.hpp"
#include <jian/nuc2d/SSTree.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(ss_tree) {
    try {
        SSTree ss_tree;
        ss_tree.make(par["seq"][0], par["ss"][0]);
        ss_tree.head->print_tree();
    } catch(const Error &e) {
        std::cout << e.what() << std::endl;
    }
}

} // namespace jian

