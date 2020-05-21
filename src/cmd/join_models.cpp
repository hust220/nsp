#include "nsp.hpp"
#include "pdb.hpp"
#include "rtsp_format.hpp"

namespace jian {

REGISTER_NSP_COMPONENT(join_models) {
    Format format;
    auto g = par.getv("global");
    int i = 0;
    Molecule mol;
    for (auto && s : g) {
        if (i >= 1) {
            auto && m = mol_read_to<Molecule>(s);
            for (auto && model : m) {
                mol.push_back(format(model));
            }
        }
        i++;
    }
    JN_OUT << mol;
}

}

