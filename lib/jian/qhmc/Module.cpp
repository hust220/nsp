#include <algorithm>
#include "Module.hpp"
#include "../utils/Factory.hpp"

namespace jian {
namespace nuc3d {
namespace quadruple {

//int Module::max_len() const {
//    return std::accumulate(d_frags.begin(), d_frags.end(), 0, [](auto &&n, auto &&frag){
//        return std::max(n, int(frag.size()));
//    });
//}

void Module::set_indices(int n, int beg, int end) {
    if (n == 0) {
        for (int i = 0; i <= end - beg; i++) {
            (*d_indices)(i, 0) = beg + i;
        }
    } else if (n == 1) {
        for (int i = 0; i <= end - beg; i++) {
            (*d_indices)(d_max_len - 1 - i, 1) = beg + i;
        }
    } else if (n == 2) {
        for (int i = 0; i <= end - beg; i++) {
            (*d_indices)(i, 2) = beg + i;
        }
    } else if (n == 3) {
        for (int i = 0; i <= end - beg; i++) {
            (*d_indices)(d_max_len - 1 - i, 3) = beg + i;
        }
    }
}

} // namespace quadruple
}
} // namespace jian


