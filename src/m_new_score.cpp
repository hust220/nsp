#include "nsp.hpp"
#include <jian/scoring/new_score.hpp>
#include <jian/utils/file.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(new_score) {
    if (par.has("l")) {
        EACH_SPLIT_LINE(par["l"][0].c_str(), " ",
            std::cout << scoring::new_score(mol_read_to<Model>(F[0])) << std::endl;
        );
    } else {
        std::cout << scoring::new_score(mol_read_to<Model>(par["s"][0])) << std::endl;
    }
}

} // namespace jian
















