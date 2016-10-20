#include <algorithm>
#include "Module.hpp"
#include "../utils/Factory.hpp"

namespace jian {
namespace qhmc {

void Module::set_indices(int n, const Frag &f, bool d = true) {
    int i, j;

    for (i = 0; i < f.size(); i++) {
        j = (d ? i : d_max_len - 1 - i);
        d_indices(j, n) = f[i];
    }
}

} // namespace qhmc
} // namespace jian


