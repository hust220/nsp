#include <algorithm>
#include "Module.hpp"
#include "../utils/Factory.hpp"

namespace jian {
namespace qhmc {

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

} // namespace qhmc
} // namespace jian


