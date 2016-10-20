#include "nsp.hpp"
#include <jian/nuc3d/ParseHelix.hpp>
#include <jian/geom.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(helix_par) {
    int n = JN_INT(par[2]);
    auto &&c = make_standard_helix(n);
    int j; for (int i = 1; i <= n; i++) {
        j = 2*n + 1 - i;
        std::cout << i << ' ' << geom::dihedral(c.row(0), c.row(i), c.row(j), c.row(2*n+1)) << std::endl;
    }
}

} // namespace jian
















