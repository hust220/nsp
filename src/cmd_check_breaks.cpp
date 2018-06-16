#include "nsp.hpp"
#include "pdb.hpp"
#include "geom.hpp"
#include "score_par_bp.hpp"

namespace jian {

static Bool are_neighbors(const Residue &r1, const Residue &r2) {
    auto && a1 = r1["O3*"];
    auto && a2 = r2["C5*"];
    Num d = geom::distance(a1, a2);
    return d < 4;
}

static Bool are_paired(const Residue &r1, const Residue &r2) {
    ParBp par(r1, r2);
    return par.is_paired();
}

REGISTER_NSP_COMPONENT(check_breaks) {
    auto g = par.getv("global");
    auto && m = mol_read_to<Model>(g[1]);
    auto rs = m.residues();
    Int l = size(rs);
    for (Int i = 0; i < l-1; i++) {
        if (!are_neighbors(rs[i], rs[i+1]) && !are_paired(rs[i], rs[i+1])) {
            JN_OUT << i+1 << ' ';
        }
    }
    JN_OUT << std::endl;
}

}

