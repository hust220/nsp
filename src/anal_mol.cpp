#include "nsp.hpp"
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

} // namespace jian

