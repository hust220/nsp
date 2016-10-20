#include <algorithm>
#include "TModule.hpp"

namespace jian {
namespace nuc3d {
namespace triple {

//int TModule::max_len() const {
//    return std::accumulate(d_frags.begin(), d_frags.end(), 0, [](auto &&n, auto &&frag){
//        return std::max(n, int(frag.size()));
//    });
//}

void TModule::set_indices(int n, int beg, int end) {
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
    }
}

} // namespace triple
}
} // namespace jian


