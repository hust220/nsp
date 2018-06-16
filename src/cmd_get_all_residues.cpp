#include <string>
#include "nsp.hpp"
#include "pdb.hpp"
#include "file.hpp"

namespace jian {

REGISTER_NSP_COMPONENT(get_all_residues) {
    auto && m = mol_read_to<Model>(par.get("s"));
    S name = file::name(par["s"][0]);
    int i = 1;
    for (auto && chain : m) {
        for (auto && res : chain) {
            mol_write(res, name + "-" + std::to_string(i) + ".pdb");
            i++;
        }
    }
}

}

