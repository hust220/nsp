#include <string>
#include "nsp.hpp"
#include "rtsp_mutate.hpp"

namespace jian {

REGISTER_NSP_COMPONENT(mutate) {
    auto global = par.getv("global");
    auto m = mol_read_to<Molecule>(global[1]);
    auto seq = global[2];
    for (auto && model : m) {
        Int i = 0;
        for (auto && chain : model) {
            for (auto && res : chain) {
                res = mutate(res, to_str(seq[i]));
                i++;
            }
        }
    }
    JN_OUT << m;
    //JN_OUT << mutate(read_model_to_chain(global[1]), global[2], "RNA");
}

}

