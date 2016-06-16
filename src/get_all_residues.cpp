#include <string>
#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/utils/file.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(get_all_residues) {
    Model m(par["pdb"][0]);
    std::string name = file::name(par["pdb"][0]);
    int i = 1;
    for (auto && chain : m) {
        for (auto && res : chain) {
            residue_to_file(res, name + "-" + std::to_string(i) + ".pdb");
            i++;
        }
    }
}

} // namespace jian

