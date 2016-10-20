#include <jian/pdb.hpp>
#include "nsp.hpp"

namespace jian {

REGISTER_NSP_COMPONENT(aa) {
    Residue r1;
    Residue r2;
    Atom a;

    a.init("P", 1, 2, 3);
    r1.push_back(a);
    r2 = r1;
}

} // namespace jian


