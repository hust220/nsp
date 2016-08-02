#include "nsp.hpp"
#include <jian/geom.hpp>
#include <jian/cg.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(score_aa) {
    auto && chain = CG1p::chain(residues_from_file(par["pdb"][0]));
    double e = 0, d;
    for (int i = 0; i < chain.size() - 1; i++) {
        d = geom::distance(chain[i][0], chain[i+1][0]);
        e += square(d - 6.1);
    }
    e /= chain.size() - 1;
    std::cout << e << ' ';

    e = 0;
    for (int i = 0; i < chain.size() - 2; i++) {
        d = geom::angle(chain[i][0], chain[i+1][0], chain[i+2][0]);
        e += square(d - 2.6);
    }
    e /= chain.size() - 2;
    std::cout << e << std::endl;
}

} // namespace jian
















