#include <algorithm>
#include "g4dna_module.hpp"
#include "factory.hpp"

BEGIN_JN
namespace qhmc {

void Module::set_indices(int n, const Frag &f, bool d = true) {
    int i, j;

    for (i = 0; i < f.size(); i++) {
        j = (d ? i : d_max_len - 1 - i);
        d_indices(j, n) = f[i];
    }
}

} // namespace qhmc
END_JN


