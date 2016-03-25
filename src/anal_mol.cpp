#include "nsp.hpp"
#include <jian/nuc3d/ParseHelix.hpp>
#include <jian/pdb/Model.hpp>
#include <jian/geom.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(anal_mol) {
    Model model(par["pdb"][0]);
    std::deque<int> ls; EACH(s, par["num"], ls.push_back(JN_INT(s)));
    std::deque<Atom> vec(ls.size());
    EACH_RES(model, auto result = std::find(ls.begin(), ls.end(), N_RES);
                    if (result != ls.end()) vec[std::distance(ls.begin(), result)] = RES["C4*"]);
    SWITCH(par[2], ("dist", SEELN(geom::distance(vec[0], vec[1]))),
                   ("ang",  SEELN(geom::angle(vec[0], vec[1], vec[2]))),
                   ("dih",  SEELN(geom::dihedral(vec[0], vec[1], vec[2], vec[3]))),
                   ("chir", SEELN(geom::chirality(vec[0], vec[1], vec[2], vec[3]))));
}

REGISTER_NSP_COMPONENT(helix_par) {
    ParseHelix parse_helix;
    int n = JN_INT(par[2]);
    auto &&c = parse_helix.make_standard_helix(n);
    int j; for (int i = 1; i <= n; i++) {
        j = 2*n + 1 - i;
        std::cout << i << ' ' << geom::dihedral(c.row(0), c.row(i), c.row(j), c.row(2*n+1)) << std::endl;
    }
}

} // namespace jian

