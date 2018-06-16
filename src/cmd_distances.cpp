#include "nsp.hpp"
#include "file.hpp"
#include "pdb.hpp"
#include "geom.hpp"

namespace jian {

static double min_distance(const Residue &r1, const Residue &r2) {
    double d, min = 99999;
    for (auto && a1 : r1) {
        for (auto && a2 : r2) {
            d = geom::distance(a1, a2);
            if (d < min) {
                min = d;
            }
        }
    }
    return min;
}

REGISTER_NSP_COMPONENT(distances) {
    auto g = par.getv("global");
    Chain c = read_model_to_chain(g[1]);
    Int l = size(c);
    for (Int i = 0; i < l; i++) {
        for (Int j = i + 1; j < l; j++) {
            JN_OUT << i+1 << ' ' << j+1 << ' ' << min_distance(c[i], c[j]) << std::endl;
        }
    }
}

}
















